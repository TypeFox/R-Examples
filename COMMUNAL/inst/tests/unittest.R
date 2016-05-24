## Contains unit tests for COMMUNAL package
test.COMMUNAL1 <- function(){
  ## Ensure proper error handling (should throw error)
  set.seed(1)
  data <- matrix(rnorm(20000, 4), ncol = 100)
  checkException(COMMUNAL(data, ks = c(1:5))) # ks must be > 1
  checkException(COMMUNAL(data, ks = c(1.5,2,4,6))) # ks must be whole numbers  
  
  data[1,1] <- -10
  checkException(COMMUNAL(data, ks = 3:5, clus.methods = c("fanny", "nmf", "sota"))) # NMF requires all entries to be positive
  
  result <- COMMUNAL(data, ks = c(2:5), seed = 1)
  checkTrue(!is.null(result$item.names))    # Check that column names are assigned if missing
  checkTrue(all(!is.na(result$measures)))  # Validation measures are computed
}

test.COMMUNAL2 <- function(){
  ## Ensure COMMUNAL runs without error in edge cases
  # One k, one clustering method, one validation measure.
  set.seed(1)
  data <- matrix(rnorm(20000, 4), ncol = 100)
  result <- COMMUNAL(data, ks = 3, clus.methods = "diana", validation = "Connectivity")
  checkEquals(dim(result$measures), c(1, 1, 1))
  
  # Two ks, one clustering method, two validation measures
  result <- COMMUNAL(data, ks = c(3,5), clus.methods = "diana", validation = c("Connectivity", "avg.silwidth"), seed = 1)
  checkEquals(length(result$ks), 2)
  checkEquals(dim(result$measures), c(2, 2, 1))
  checkEquals(dim(result$getClustering(3)), c(100, 1))
  checkException(result$getClustering(4)) # Invalid k.
}