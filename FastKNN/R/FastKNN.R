#'k-Nearest Neighbors
#'the \code{k.nearest.neigbors} gives the list of points (k-Neigbours) that are closest
#'to the row i in descending order.
#'@param i is from the numeric class and is a row from the distance_matrix.
#'@param distance_matrix is a nxn matrix.
#'@param k is from the numeric class and represent the number of neigbours that the function will return.
#'@return a k vector with the k closest neigbours to the i observation.
#'@details
#'The output of this function is used in the \code{knn_test_function} function.
#'@seealso \code{order}
#'@export
#'@import assertthat

k.nearest.neighbors <- function(i, distance_matrix, k = 5)
{
  #Test the inputs
  not_empty(i); not_empty(distance_matrix);not_empty(k);
  
  ordered_neighbors <- order(distance_matrix[i, ]) 
  return(ordered_neighbors[2:(k + 1)])
}

#'KNN Test
#'The knn_test_function returns the labels for a test set using the k-Nearest Neighbors Clasification method.
#'@param dataset is a matrix with the features of the training set
#'@param test is a matrix where the columns are the features of the test set
#'@param distance is a nxn matrix with the distance between each observation of the test set
#'and the training set
#'@param labels is a nx1 vector with the labels of the training set
#'@param k is from the numeric class and represent the number of neigbours to be use in the classifier.
#'@return a k vector with the predicted labels for the test set.
#'@seealso \code{k.nearest.neighbors}
#'@seealso \code{sample}
#'@export
#'@import assertthat
#'@examples
#' # Create Data for restaurant reviews
#' training <- matrix(rexp(600,1), ncol=2)
#' test  <- matrix(rexp(200,1), ncol=2)
#' # Label "Good", "Bad", "Average"
#' labelsExample <- c(rep("Good",100), rep("Bad",100), rep("Average",100))
#' # Distance Matrix
#' distanceExample<-Distance_for_KNN_test(test, training)
#' # KNN
#' knn_test_function(training, test, distanceExample,labelsExample, k = 3)

knn_test_function <- function(dataset, test, distance,labels, k = 3)
{
  #Test the inputs
  not_empty(dataset); not_empty(test);not_empty(distance);
  not_empty(labels); not_empty(k);
  if(ncol(dataset)!=ncol(test))
    stop("The test set and the training set have a different number of features")
  
  predictions <- rep(0, nrow(test))
 for (i in 1:nrow(test)){
   indices <- k.nearest.neighbors(i, distance, k = k)
   vec_mas_cer <- table(labels[indices])
   predictions[i] <- sample(names(which(vec_mas_cer == max(vec_mas_cer))),size = 1)}
 return(predictions)}


#'KNN Training
#'The knn_training_function returns the labels for a training set using the 
#'k-Nearest Neighbors Clasification method.
#'@param dataset is a matrix with the features of the training set
#'@param distance is a nxn matrix with the distance between each observation of the training set
#'@param label is a nx1 vector with the labels of the training set
#'@param k is from the numeric class and represent the number of neigbours to be use in the classifier.
#'@return a k vector with the predicted labels for the training set.
#'#'@details
#'This function is use to check the quality of the Classifier. Because then the predicted labels
#'are compared with the true labels
#'@seealso \code{k.nearest.neighbors}
#'@seealso \code{sample}
#'@export
#'@import assertthat

knn_training_function <- function(dataset, distance,label, k = 1)
{
  #Test the inputs
  not_empty(dataset); not_empty(distance);
  not_empty(labels); not_empty(k);
  if(isSymmetric(distance)==FALSE)
    stop("The distance Matrix is not symmmetric")
  
predictions <- rep(0, nrow(dataset))
 # For every point in the dataset
 for (i in 1:nrow(dataset)){
   indices <- k.nearest.neighbors(i, distance, k = k)
   vec_mas_cer <- table(label[indices])
   predictions[i] <- sample(names(which(vec_mas_cer == max(vec_mas_cer))), size = 1)}
 return(predictions)}

#'Distance for KNN Test
#'The Distance_for_KNN_test returns the distance matrix between the test set and the training set.
#'@param train_set is a matrix with the features of the training set
#'@param test_set is a matrix where the columns are the features of the test set
#'@return a distance matrix
#'@seealso \code{knn_test_function}
#'@seealso \code{pdist}
#'@export
#'@import pdist
#'@import assertthat


Distance_for_KNN_test<-function(test_set, train_set)
{
  #Test the inputs
  not_empty(test_set); not_empty(train_set);
  if(ncol(test_set)!=ncol(train_set))
    stop("The test set and the training set have a different number of features")
  
  distanciapost<-matrix(NA,nrow(test_set),nrow(train_set))
  
  for (i in 1:nrow(test_set)){
    distanciapost[i,]<-as.matrix(pdist(train_set,test_set[i,]))}
  return(distanciapost)
}

