summary.transformedTrain <- function(object,...) {

  Y <- seq(object)
  nbSpikes <- length(object)
  M <- Y - as.numeric(object)
  uniformOnTTime <- c(in95=all(-1.36*sqrt(nbSpikes) <= M & M <= 1.36*sqrt(nbSpikes)),
                      in99=all(-1.63*sqrt(nbSpikes) <= M & M <= 1.63*sqrt(nbSpikes))
                      )

  lambda <- 1-exp(-diff(object))
  M <- (1:(nbSpikes-1))/(nbSpikes-1)-sort(lambda)
  BermanTest <- c(in95=all(-1.36/sqrt(nbSpikes-1) <= M & M <= 1.36/sqrt(nbSpikes-1)),
                  in99=all(-1.63/sqrt(nbSpikes-1) <= M & M <= 1.63/sqrt(nbSpikes-1))
                  )
  vt <- varianceTime(object,...)
  vtResult <- sapply(seq(vt$CI),
                     function(idx) {
                       nb <- length(vt$s2)
                       outNB <- nb-sum(vt$ciLow[idx,] <= vt$s2 &
                                       vt$s2 <= vt$ciUp[idx,])
                       pbinom(outNB,nb,1-vt$CI[idx])
                     }
                     )
  names(vtResult) <- paste(vt$CI)

  list(uniformOnTTime=uniformOnTTime,
       BermanTest=BermanTest,
       VarTime=vtResult)
  
}
