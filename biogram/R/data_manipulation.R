# #' Prepare data for further analysis
# #'
# #' Creates random stratified subsample of population.
# #' 
# #' @param data_set data set.
# #' @param pos indices of positive data.
# #' @param neg indices of negative data.
# #' @param train_size a vector of length 2 - size of positive and negative training
# #' subsample.
# #' @param test_size a vector of length 2 - size of positive and negative test
# #' subsample.
# #' @return a list of a length 2 containing cases belonging to respectively train and test 
# #' sample.
# #' @details Check it for 0 sample size.
# #' @export
# #' @examples 
# #' #data(iris)
# #' #prepare_data(iris, 1:20, 101:120, c(5, 5), c(8, 8))
# 
# prepare_data <- function(data_set,
#                          pos, 
#                          neg, 
#                          train_size, #n pos, n neg
#                          test_size) {
#   id_pos <- data_sample(pos, train_size[1], test_size[1])
#   id_neg <- data_sample(neg, train_size[2], test_size[2])
#   train_dat <- data_set[c(id_pos[["train"]], id_neg[["train"]]), ]
#   test_dat <- data_set[c(id_pos[["test"]], id_neg[["test"]]), ]
#   list(train = train_dat,
#        test = test_dat)
# }
# 
# 
# #' Subsample data
# #'
# #' Creates random subsample of population using indices.
# #'
# #' @param indices of data.
# #' @param train_size size of training sample.
# #' @param test_size size of test sample.
# #' @return a list of a length 2 containing train and test sample indices.
# #' @details Check it for 0 sample size.
# #' @export
# #' @examples data_sample(1:20, 5, 5)
# 
# data_sample <- function(indices, train_size, test_size) {
#   #sampled indices  
#   s_indices <- sample(indices, train_size + test_size)
#   train_ind <- s_indices[0L:train_size]
#   if (train_size != 0) {
#     test_ind <- s_indices[-c(0L:train_size)][0L:test_size]
#   } else {
#     test_ind <- s_indices[0L:test_size]
#   }
#   list(train = train_ind, test = test_ind)
# }
