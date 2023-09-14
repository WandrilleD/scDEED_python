#!/usr/bin/python3

## this code proposes a python implementation of the method originally described in:
##
## "scDEED: a statistical method for detecting dubious 2D single-cell embeddings "
## Lucy Xia, Christy Lee, and Jingyi Jessica Li
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10168265/
##
## the original R implementation can be found at https://github.com/JSB-UCLA/scDEED/



import numba
from umap import UMAP
from sklearn.manifold import TSNE
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
import scipy.stats as stats


def get_reliability_score(embed ,
                          high_neighbors_indexes,
                          fraction_selected = 0.5):
    """
    Computes the reliability score of each points in the embedded dataset.
    The reliability score is, for each point, the Pearson correlation of 2 distance vectors:
     * euclidean distance in the low-dimensional space to the x% closest neighbors in low-dimensional space
     * euclidean distance in the low-dimensional space to the x% closest neighbors in high-dimensional space
    
    further definitions:
     * n : number of points
     * m : number of dimensions in the low-dimensional space
    
    Takes:
     - embed (n*m numpy.array) : low-dimensional representation of the dataset 
     - high_neighbors_indexes (n*n numpy.array) : matrix where row i corresponds to indexes of point i neighbors sorted by their euvclidean distance in the high dimensional space
     - fraction_selected (float , default = 0.5) : fraction of closest neighbors to consider to compute the reliability score
    
    Returns
     (numpy.array) : reliabilty score of each point
    """

    ## order neighbors in the low dim space
    dist_low = squareform( pdist(embed) )
    low_neighbors_indexes = list( map( np.argsort , dist_low ) )

    numberselected = int(embed.shape[0]*fraction_selected)

    r_embed = []

    for i in range(embed.shape[0]):
        r_embed.append( stats.pearsonr( dist_low[i, high_neighbors_indexes[i][ 1:(numberselected+1) ] ],
                                       dist_low[i, low_neighbors_indexes[i][ 1:(numberselected+1) ] ]).statistic )
    return np.array(r_embed)




def get_norm_score( test, ref ):
    """
    Takes:
     - (numpy.array) score for a test dataset
     - (numpy.array) score for a reference (ie. null) dataset
     
    Returns
        (numpy.array) : for each element in the test dataset: which fraction of the reference set have a lower score
    """
    O = np.argsort( np.concatenate((test,ref) ) )
    N = len(ref)
    scores = [0.0]*N
    current_score = 0.0
    score_steps = 1/N
    for o in O:
        if o<N:
            scores[o] = current_score
        else:
            current_score += score_steps

    return np.array(scores)

def norm_score_to_scDEED(score):
    assign = np.array(['other']*score.shape[0] , dtype='<U7')
    assign[score < 0.05] = 'dubious'
    assign[score > 0.95] = 'good'
    return assign

    

def get_neighbors_matrix(data):
    dist_data = squareform( pdist(data) )
    return np.array( list( map( np.argsort , dist_data ) ) )


def permute_data(data):
    """row permutation for each column"""
    P = np.zeros(arr.shape)

    for i in range( arr.shape[1] ):
        P[:,i] = arr[rng.permutation(arr.shape[0]),i]

    return P
    

def scDEED( data , 
           embedding_function, 
           data_permuted = None,
           high_neighbors_indexes = None,
           high_neighbors_permuted_indexes = None,
           fraction_selected = 0.5):
    
    
    if data_permuted is None:
        ## create a permuted version of the data
        data_permuted = permute_data(data)
        ## order neighbors in the permuted high dim space
        high_neighbors_permuted_indexes = get_neighbors_matrix(data_permuted)
        
    if high_neighbors_permuted_indexes is None:
        ## order neighbors in the permuted high dim space
        high_neighbors_permuted_indexes = get_neighbors_matrix(data_permuted)
    
    if high_neighbors_indexes is None:
        ## order neighbors in the original high dim space
        high_neighbors_indexes = get_neighbors_matrix(data)
    
    
    ## embed orginal data
    embed = embedding_function( data )
    
    ## embed permuted data
    embed_perm = embedding_function( data_permuted )
    
    r_embed = get_reliability_score(embed , high_neighbors_indexes, fraction_selected)
    r_embed_permuted = get_reliability_score(embed_perm , high_neighbors_permuted_indexes, fraction_selected)
    
    norm_score = get_norm_score( r_embed, r_embed_permuted )
    
    scDEED_score = norm_score_to_scDEED(norm_score)
    
    return {"embedded_data":embed , 'normalized_score':norm_score , 'scDEED_score':scDEED_score}
