optimal_coding <-
function (data, dimens, alpha) {

   nCandidates <- dim(data)[2] - 1
   tData <- t(data)

   bestCandidate <- 1
   bestLogMargLik <- findLogMargLik (c(1,nCandidates + 1), tData, dimens, alpha)
   
   for (i in 2:nCandidates) {
     logMargLik <- findLogMargLik (c(i,nCandidates + 1), tData, dimens, alpha)
     if (logMargLik > bestLogMargLik) {
       bestCandidate <- i
       bestLogMargLik <- logMargLik
     }
   }

   return(list(data[,bestCandidate], dimens[bestCandidate]))
}
