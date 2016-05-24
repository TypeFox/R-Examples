# test for equality of the saved result and the actual result of a dontrun example
test.qqplot1 <- function() {
  # we compute the actual result
  P <- Pois(5)
  B <- Binom(size=2000,prob=5/2000)
  runit.dontrunQQPlot1.actual <- qqplot(B,P, nosym.pCI=TRUE)
  
  # we load the stored result
  #   we assume that this test is called from within the script in the upper directory
  load("unitTests/runit.dontrunQQPlot1.save")
  
  # we compare the stored result with the calculated one
  #   (a comparison with identical (ignoring the environment) gives FALSE...
  result <- all.equal(runit.dontrunQQPlot1.actual,
                      runit.dontrunQQPlot1.save)
  
  # we check whether the result is TRUE and if not, we write the message
  #   coming from the result
  checkEquals(is.logical(result) && result, TRUE, msg=paste(result, sep="", collapse="\n"))
}


# test for equality of the saved result and the actual result of a dontrun example
test.qqplot2 <- function() {
  # we compute the actual result
  mylist <- UnivarLebDecDistribution(discretePart=Binom(3,.3), acPart=Norm(2,2),
                                     acWeight=11/20)
  mylist2 <- mylist+0.1
  
  runit.dontrunQQPlot2.actual <- qqplot(mylist,mylist2,nosym.pCI=TRUE)
  
  # we load the stored result
  #   we assume that this test is called from within the script in the upper directory
  load("unitTests/runit.dontrunQQPlot2.save")
  
  # we compare the stored result with the calculated one
  #   (a comparison with identical (ignoring the environment) gives FALSE...
  result <- all.equal(runit.dontrunQQPlot2.actual,
                      runit.dontrunQQPlot2.save)
  
  # we check whether the result is TRUE and if not, we write the message
  #   coming from the result
  checkEquals(is.logical(result) && result, TRUE, msg=paste(result, sep="", collapse="\n"))
}


# test for equality of the saved result and the actual result of a dontrun example
test.qqplot3 <- function() {
  # we compute the actual result
  mylist3 <- UnivarMixingDistribution(Unif(0,0.3),Unif(0.6,1),mixCoeff=c(0.8,0.2))
  mylist4 <- UnivarMixingDistribution(Unif(0,0.3),Unif(0.6,1),mixCoeff=c(0.6,0.4))
  
  runit.dontrunQQPlot3.actual <- qqplot(mylist3,mylist4,nosym.pCI=TRUE)
  
  # we load the stored result
  #   we assume that this test is called from within the script in the upper directory
  load("unitTests/runit.dontrunQQPlot3.save")
  
  # we compare the stored result with the calculated one
  #   (a comparison with identical (ignoring the environment) gives FALSE...
  result <- all.equal(runit.dontrunQQPlot3.actual,
                      runit.dontrunQQPlot3.save)
  
  # we check whether the result is TRUE and if not, we write the message
  #   coming from the result
  checkEquals(is.logical(result) && result, TRUE, msg=paste(result, sep="", collapse="\n"))
}

