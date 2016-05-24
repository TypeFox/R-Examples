## the body of this function is the example for single values in data.tex

testSd <- function(sdev, eps ) {
    Wt <- Inf
    testSd <- try(min(sdev) > eps)
    if(identical(testSd, TRUE))
      Wt <- 1/sdev
    else if(!identical(testSd, FALSE)) {
        if(is(testSd, "try-error"))
          stop("Encountered error in testing sdev \"",
               testSd, "\"")
        else
          stop("Testing sdev produced an invalid result: ",
               summaryString(testSd))
    }
    Wt
}
