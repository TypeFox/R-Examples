# function that generates a by-cluster correlation matrix.

clusterCorr <- function(observed_cor_matrix,cluster_vector) {
	num_vertices = nrow(observed_cor_matrix)

	cluster_cor_mat <- observed_cor_matrix
	for (i in 1:num_vertices) {
		for (j in 1:num_vertices) {
		cluster_cor_mat[i,j] = mean(observed_cor_matrix[which(
		  cluster_vector[row(observed_cor_matrix)]==cluster_vector[i] & 
		  cluster_vector[col(observed_cor_matrix)]==cluster_vector[j])])
		}
	}
	
	return(cluster_cor_mat)
}

clustConfigurations <- function(vertices,hclustresult,observedcorrelation)
{

	resultlist <- list()
	correlations <- vector()
	for (i in 1:(vertices)) {
		cluster_result <-list(label=NA,clusters=NA,correlation=NA)
		# Generate the clusters for this solution.
		cluster_result$label <- paste('number of clusters: ', i)
		clusters <- cutree(hclustresult, k = i)
		cluster_result$clusters <- clusters
	
		# Compute within- and between-cluster correlations, using
		# the function we defined above.
		cluster_cor_mat <- clusterCorr(observedcorrelation, clusters)

		# Correlate the by-cluster correlation matrix with the 
		# observed correlation matrix using gcor() from the 
		# sna package.
		clustered_observed_cors <- gcor(cluster_cor_mat, observedcorrelation)
	  
		cluster_result$correlation <-(clustered_observed_cors)
		resultlist <-c(resultlist,cluster_result)
		correlations<-c(correlations,clustered_observed_cors)
	}
		#returnlist<-list(output=resultlist,correlations=correlations)
		plot(correlations, type = "b")
		resultlist$correlations <- correlations
	return(resultlist)
}

permute_matrix <- function(mem_vector, adj_matrix){
	adj_matrix_with_comms <- cbind(mem_vector, adj_matrix)

	# Sort rows by community.
	adj_matrix_sorted_rows <- adj_matrix_with_comms[order(adj_matrix_with_comms[,1],decreasing=TRUE),]

	# Remove the column listing the communities to get us back to
	# a square matrix.
	adj_matrix_sorted_rows <- adj_matrix_sorted_rows[,-1]
	
	# Transpose the matrix so it's now the columns that are
	# sorted.
	adj_matrix_sorted_cols <- t(adj_matrix_sorted_rows)
	
	# Run through the above steps one more time to sort the 
	# (now unsorted) rows. 
	adj_matrix_with_comms <- cbind(mem_vector, adj_matrix_sorted_cols)
	adj_matrix_sorted_rows <- adj_matrix_with_comms[order(adj_matrix_with_comms[,1],decreasing=TRUE),]
	adj_matrix_sorted <- adj_matrix_sorted_rows[,-1]
	return(adj_matrix_sorted)

}


generate_cluster_cor_mat <- function(observed_cor_matrix,
		cluster_vector) {
	num_vertices = nrow(observed_cor_matrix)
	
	cluster_cor_mat <- observed_cor_matrix
	for (m in 1:num_vertices) {
		for (n in 1:num_vertices) {
			cluster_cor_mat[m,n] = mean(observed_cor_matrix[which(
									cluster_vector[row(observed_cor_matrix)]==cluster_vector[m] & 
											cluster_vector[col(observed_cor_matrix)]==cluster_vector[n])])
		}
	}
	return(cluster_cor_mat)
}
