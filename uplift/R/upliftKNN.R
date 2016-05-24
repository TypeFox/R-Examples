######################################################################
# Uplift K-Nearest Neighbour Regression and Classification
######################################################################

# train = matrix or data frame of training set cases.
# test = matrix or data frame of test set cases. A vector will be interpreted as a row vector for a single case.
# y = a numeric response (must be coded as 0/1 for binary response)
# ct = factor or numeric vector representing the treatment to which each train case is assigned. 
#      At least 2 groups is required (e.g. treatment and control). Multi-treatments is
#      also supported
# k = number of neighbours considered.
# dist.method = the distance to be used in calculating the neighbors. Any method supported in function
#               \link{dist} is valid.
# p = the power of the Minkowski distance.
# ties.meth = method to handle ties for the kth neighbor. The default is "min" which uses all
#             ties. Alternatives include "max" which uses none if there are ties for the k-th
#             nearest neighbor, "random" which selects among the ties randomly and "first"
#             which uses the ties in their order in the data.
# agg.method = method to combine responses of the nearest neighbors, defaults to "majority"
#              for classification and "mean" for continuous responses
#              majority: with ties broken at random
# In details, the logic for the code follows closely that the knn and xknnflex package, the later 
   #currently discontinued in CRAN.

### majority vote function

majority <- function(x) {
  x <- as.factor(x)
  n <- nlevels(x)
  votes <- rep(0, n)
  for (i in 1:length(x)) votes[as.integer(x[i])] <- votes[as.integer(x[i])] + 1
  as.numeric(levels(x)[order(votes, decreasing=TRUE, sample(1:n,n))[1]])
}

### uplift k-nearest neighbour

upliftKNN <- function(train, test, y, ct, k = 1, dist.method = "euclidean", p = 2, 
                      ties.meth = "min", agg.method = "mean") {
    
  ### perform various checks
  train <- as.matrix(train)
  if(is.null(dim(test))) dim(test) <- c(1, length(test))
  test <- as.matrix(test)
  if(any(is.na(train)) || any(is.na(test)) || any(is.na(y)) || any(is.na(ct)))
    stop("upliftKNN: no missing values are allowed")
  ptr <- ncol(train)
  ntr <- nrow(train)
  fct <- as.factor(ct)
  if (length(unique(fct)) < 2) stop("uplift: ct must have at least 2 distinct values")
  if (!is.numeric(y)) stop("uplift: y must be a numeric variable")
  if (length(y) != ntr) stop("uplift: 'train' and 'y' have different lengths")
  if (length(ct) != ntr) stop("uplift: 'train' and 'ct' have different lengths")
  if(ntr < k) {
    warning(gettextf("uplift: k = %d exceeds number %d of patterns. Reset k = %d.", k, ntr, ntr),
            domain = NA)
    k <- ntr
  }
  if (k < 1)
    stop(gettextf("uplift: k = %d must be at least 1", k), domain = NA)
  nte <- nrow(test)
  if(ncol(test) != ptr) stop("uplift: dims of 'test' and 'train' differ")
  am <- charmatch(tolower(agg.method), c("mean", "majority"))
  if (is.na(am)) stop("uplift: agg.method must be one of 'mean' or 'majority'")
      
  ### compute distance matrix
  x <- rbind(train, test)
  dist. <- as.matrix(dist(x, dist.method, p))
  
  ### only need the rows for the train data and columns for the test data
  dist. <-  as.data.frame(dist.[1:ntr, ntr + 1:nte])
  
  ### split distance matrix by the number of treatment grups
  dist.split <- split(dist., fct, drop = TRUE)
  y.split <- split(y, ct, drop = TRUE)
  ranks <- lapply(1:length(dist.split), function(i) t(apply(dist.split[[i]], 2, 
                  function(x) rank(x, ties.method = ties.meth))))
  agg <- lapply(1:length(ranks), function(i)
               apply(ranks[[i]], 1, function(x) 
                 apply(data.frame(y.split[[i]][x <= k]), 2, agg.method)))  
  
  # create output matrix with outcomes associated with each test obs and treatment type
  res <- matrix(nrow = nte, ncol = length(unique(ct)))
  for (i in 1:length(unique(ct))) {
    res[, i] <- agg[[i]]
  }
  
  colnames(res) <- names(dist.split)
  return(res)
  
}



