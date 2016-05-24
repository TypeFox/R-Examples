predictLogProb <- function (object,
                            newdata
                            ) {

  if (!inherits(object, c("ssanova", "ssanova0"))) 
    stop("object should be a ssanova or a ssanova0 obbject.")
  allTerms <- object$terms$labels[-1]
  isInteraction <- sapply(strsplit(allTerms,":"),length) > 1
  dfTerms <- allTerms[!isInteraction]
  dataRange <- sapply(dfTerms, function(n) range(object$mf[,n]))
  for (n in dfTerms) {
    tooSmall <- newdata[,n] < dataRange[1,n]
    newdata[tooSmall,n] <- dataRange[1,n]
    tooBig <- newdata[,n] > dataRange[2,n]
    newdata[tooBig,n] <- dataRange[2,n]
  } ## end of loop on n
  newPred <- predict(object, newdata = newdata)
  sum(with(newdata, event * newPred) - log(1 + exp(newPred)))

}
