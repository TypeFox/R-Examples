`featureSelect` <-
  function(X, y, numFeatures=10){
    stopifnot(is.matrix(X), is.numeric(X), is.factor(y))
    if(nrow(X)!=length(y)) stop("X must be a data matrix with number of
                                rows equal to the length of y")
    ans <- fitted(aov(X ~ y))
    dim(ans) <- dim(X)
    SSB <- apply(ans, 2, function(x) sum((x - mean(x))^2))
    order(SSB, decreasing=TRUE)[1:numFeatures]
  }
