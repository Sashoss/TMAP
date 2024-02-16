
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics.cluster import calinski_harabasz_score


class CLUSTER:
    def __init__(self, components, method):
        

    def cluster_genes_Kmeans(components, n_clusters_range):
        """Clusters the genes in the given DataFrame using the given factor columns and returns the most optimal number of clusters.

        Args:
            components:  Pandas DataFrame with columns containing the factors/components.
            n_clusters_range: A range of integers representing the number of clusters to test.

        Returns:
            The most optimal number of clusters, as determined by the silhouette score.
        """
        # Extract the components as a 2D numpy array
        components_array = components.values

        # Evaluate the clustering results for each number of clusters
        silhouette_scores = []
        cluster_dict = dict()
        for n_clusters in n_clusters_range:
            # Cluster the genes using the k-means algorithm
            kmeans = KMeans(n_clusters=n_clusters)
            clusters = kmeans.fit_predict(components_array)

            # Compute the silhouette score for the clusters
            s_score = silhouette_score(components_array, clusters)
            silhouette_scores.append(s_score)
            cluster_dict[s_score] = clusters

        # Find the number of clusters with the highest silhouette score
        optimal_n_clusters = n_clusters_range[np.argmax(silhouette_scores)]

        return optimal_n_clusters, cluster_dict[optimal_n_clusters]


    def heirarchical_clustering(components, n_clusters_range):

        # Perform hierarchical clustering on the components
        components_array = components.values
        Z = linkage(components_array, method='ward')
        # Determine the optimal number of clusters using the Calinski-Harabasz score
        scores = []
        for num_clusters in range(2, n_clusters_range):
            clusters = fcluster(Z, num_clusters, criterion='maxclust')
            score = calinski_harabasz_score(components_array, clusters)
            scores.append(score)

        optimal_num_clusters = np.argmax(scores) + 2

        # Use the optimal number of clusters to assign samples to clusters
        clusters = fcluster(Z, optimal_num_clusters, criterion='maxclust')

        return clusters
    
    def check_cluster_similarity(self, clustering_1, clustering_2):
        # Check that the clustering done here is similar to user asigned sample annotation
        if len(clustering_1) != len(clustering_2):
            return False
  
        # Create a mapping from clusters in clustering_1 to clusters in clustering_2
        cluster_mapping = {}
        for i in range(len(clustering_1)):
            cluster_1 = clustering_1[i]
            cluster_2 = clustering_2[i]
            if cluster_1 not in cluster_mapping:
                cluster_mapping[cluster_1] = cluster_2
            elif cluster_mapping[cluster_1] != cluster_2:
                # If a cluster in clustering_1 is mapped to multiple clusters in clustering_2, the clusterings are not similar
                return False
  
        # If all checks pass, the clusterings are similar
        return True


# add a function that takes sample ids and returns the clusters and gives some idea of how cluster looks as compared to the groups assigned
# take this heirarchical clusters do stuff -> group specific markers
