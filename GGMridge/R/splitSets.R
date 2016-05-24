splitSets <- function(cv, i, x) {

  testindex <- cv[cv[,1L] == i,2L]
  lt <- length(testindex)
  if( lt < 0.5 ) {
    stop("Could not divide data according to specified fold.", 
    call. = FALSE)
  }

  train <- x[-testindex,,drop = FALSE]
  train <- scale(x = train, center = TRUE, scale = TRUE)

  test <- x[testindex,,drop = FALSE]
  tst <- lt > 1.5
  test <- scale(x = test, center = tst, scale = tst)

  return(list("train" = train,
              "test" = test))

}
