
size_of_rows <- 5
size_of_columns <- 10
n_clusters <- 3
my_medoids <- NULL

# test w ff input data
test.correct_ff_mediod_names <- function(){
	
# Generate random data
input_dataset <- matrix(rnorm(size_of_rows * size_of_columns), ncol=size_of_columns)
rnames <- c('v','w','x','y','z')
rownames(input_dataset) <- rnames
	
# Create different object types accepted by parallel pam
 data_symmetric_matrix <- 1-cor(t(input_dataset))
   data_distance_matrix <- as.dist(data_symmetric_matrix)
	data_binary_file <- ff(data_symmetric_matrix, vmode="double", dim=c(size_of_rows, size_of_rows),
						   dimnames=list(rnames,rnames))
    
# Execute original version
   original_pam_result <- pam(data_distance_matrix, n_clusters)
	
# Execute parallel version passing different object types on input
   ppam_result_binary_file <- ppam(data_binary_file, n_clusters)
	
# Compare medoids
	
   invisible(checkEquals(original_pam_result$medoids,
                                ppam_result_binary_file$medoids, " Medoids should have same names."))
	
}

test.symmetric_correct_mediod_names <- function(){
	
# Generate random data
input_dataset <- matrix(rnorm(size_of_rows * size_of_columns), ncol=size_of_columns)
rownames(input_dataset) <- c('v','w','x','y','z')

# Create different object types accepted by parallel pam
data_symmetric_matrix <- 1-cor(t(input_dataset))
data_distance_matrix <- as.dist(data_symmetric_matrix)

# Execute original version
original_pam_result <- pam(data_distance_matrix, n_clusters)

# Execute parallel version passing different object types on input
ppam_result_symmetric_matrix <- ppam(data_symmetric_matrix, n_clusters)

# Compare medoids

invisible(checkEqualsNumeric(original_pam_result$medoids,
ppam_result_symmetric_matrix$medoids, "Symmetric medoid labels are the same."))

invisible(checkEqualsNumeric(original_pam_result$id.med,
ppam_result_symmetric_matrix$id.med, "Symmetric meds are the same."))
}



test.distance_correct_mediod_names <- function(){
	
# Generate random data
    input_dataset <- matrix(rnorm(size_of_rows * size_of_columns), ncol=size_of_columns)
	rownames(input_dataset) <- c('v','w','x','y','z')
# Create different object types accepted by parallel pam
    data_symmetric_matrix <- 1-cor(t(input_dataset))
    data_distance_matrix <- as.dist(data_symmetric_matrix)
# Execute original version
    original_pam_result <- pam(data_distance_matrix, n_clusters)
	
# Execute parallel version passing different object types on input
    ppam_result_distance_matrix <- ppam(data_distance_matrix, n_clusters)
	
# Compare medoids
    invisible(checkEqualsNumeric(original_pam_result$medoids,
                                 ppam_result_distance_matrix$medoids, "Distance medoid labels are the same."))

	invisible(checkEqualsNumeric(original_pam_result$id.med,
		ppam_result_distance_matrix$id.med, "Distance med labels are the same."))
	
	}