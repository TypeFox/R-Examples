plotcvriskMSM <- function(cvriskMSMobject, type = "all"){
  cvm <- attributes(cvriskMSMobject)$cvpl.matrix
  mcvm <- apply(cvm, MARGIN = 1, FUN = mean)
  if(type == "all"){
    plot(cvm[, 1], type="n", xlab="Number of boosting iterations", 
         ylab="Cox Partial Likelihood", ylim = range(cvm))
    apply(cvm, MARGIN = 2, FUN = lines, col = "lightgrey")
    lines(mcvm)
    lines(rep(which.max(mcvm), 2), c(-1e5, mcvm[which.max(mcvm)]), lty=2)
  }else{
    plot(mcvm, type="l", xlab="Number of boosting iterations", 
         ylab="Mean(Cox Partial Likelihood)")
    lines(rep(which.max(mcvm), 2), c(-1e5, mcvm[which.max(mcvm)]), lty=2)
    }
}