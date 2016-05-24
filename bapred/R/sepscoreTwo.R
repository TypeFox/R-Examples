sepscoreTwo <-
function(xb1, xb2, k=10) {

  if(k >= (nrow(xb1) + nrow(xb2)))
    k <- (nrow(xb1) + nrow(xb2))-1

  Xall <- rbind(xb1, xb2)
  batchonetwo <- c(rep(1, nrow(xb1)), rep(2, nrow(xb2)))
  
  firstbatch <- names(table(batchonetwo))[1]
  secondbatch <- names(table(batchonetwo))[2]
  
  Xfirst <- Xall[batchonetwo == firstbatch,]
  Xsecond <- Xall[batchonetwo == secondbatch,]
  
  distances <- mapply(function(x, y) sqrt(rowSums(sweep(Xall[-y,], 2, x, "-")^2)), 
    data.frame(t(Xfirst)), which(batchonetwo == firstbatch))
  nearneigh <- mapply(function(x, y) sum(batchonetwo[-y][order(x)[1:k]]!=firstbatch), 
    data.frame(distances), which(batchonetwo == firstbatch))

  mixscorefirst <- abs(sum(nearneigh)/(k*nrow(Xfirst)) - sum(batchonetwo != firstbatch)/(nrow(Xall)-1))
  
  
  distances <- mapply(function(x, y) sqrt(rowSums(sweep(Xall[-y,], 2, x, "-")^2)), 
    data.frame(t(Xsecond)), which(batchonetwo == secondbatch))
  nearneigh <- mapply(function(x, y) sum(batchonetwo[-y][order(x)[1:k]]!=secondbatch), 
    data.frame(distances), which(batchonetwo == secondbatch))

  mixscoresecond <- abs(sum(nearneigh)/(k*nrow(Xsecond)) - sum(batchonetwo != secondbatch)/(nrow(Xall)-1))

  mean(c(mixscorefirst, mixscoresecond))  

}
