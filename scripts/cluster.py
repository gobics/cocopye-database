import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering
import pandas as pd
from itertools import compress


def agglo_half_manhatten_clustering(count_mat, distance_threshold):
    transformed_count_mat = 0.5 * (count_mat.T / count_mat.sum(axis=1)).T
    clust = AgglomerativeClustering(
        n_clusters=None, distance_threshold=distance_threshold, linkage="complete", metric="manhattan")
    labs_agg = clust.fit_predict(transformed_count_mat)
    return labs_agg, transformed_count_mat


def cluster(dataframe):
    count_mat = dataframe.to_numpy()
    index = dataframe.index.to_list()

    distance_threshold = 0.1

    cluster_labels, transformed_count_mat = agglo_half_manhatten_clustering(
                count_mat, distance_threshold)

    ulabs, uinvs, ucnts = np.unique(
        cluster_labels, return_inverse=True, return_counts=True)

    new_index = list(compress(index, (ucnts == 1)[uinvs]))
    x_ref1_mat = count_mat[(ucnts == 1)[uinvs]]
    x_ref2_mat = np.zeros(((ucnts > 1).sum(), transformed_count_mat.shape[1]))

    for i, lab in enumerate(ulabs[ucnts > 1]):
        sum_vec = np.sum(pairwise_distances(
            transformed_count_mat[cluster_labels == lab], metric='manhattan'), axis=0)
        x_ref2_mat[i] = count_mat[cluster_labels == lab][sum_vec.argmin()]
        new_index.append(list(compress(index, cluster_labels == lab))[
                         sum_vec.argmin()])

    x_ref_mat = np.concatenate((x_ref1_mat, x_ref2_mat), axis=0)

    columns = []
    for i in range(1, 17127):
        columns.append(f"PF{'0' * (5-len(str(i)))}{i}")
    df = pd.DataFrame(x_ref_mat, index=new_index,
                      dtype="int", columns=columns)
    return df
