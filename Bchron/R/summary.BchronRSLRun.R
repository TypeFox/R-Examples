summary.BchronRSLRun <-
function(object,...) {
  cat('95% posterior intervals...\n')
  cat('Power Lower Upper \n')
  pow.names = c('mean','linear','quadratic','cubic','quartic','quintic')  
  for(j in 1:(object$degree+1)) {
    cat(pow.names[j],round(quantile(object$samples[,j],probs=c(0.025,0.975)),4),'\n')
  } 
}
