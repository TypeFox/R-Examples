#==============================================================================
#' @importFrom graphics plot hist box axis mtext curve
#' @importFrom stats dnorm 
MixModelPlot = function (model, main1="", xlab1="Iteration", 
                         ylab1="Log-likelihood", col1="darkblue", main2="", 
                         xlab2="Data", ylab2="Density", breaks="FD", 
                         col2=NULL, loglik=TRUE, labels=TRUE, ...) 
{
  # Plot the log-likelihoods and the estimated 
  # univerate Gaussian/Student mixture model
  # 
  # Date: 
  #   Revised: February 15, 2015
  #
  # Author: Jiarui Ding <jiaruid@cs.ubc.ca>
  #   Department of Computer Science, UBC
  #   Department of Molecular Oncology, BC Cancer Agency 
  #
    
  if(loglik == TRUE) {
    plot(model$loglik, main=main1, xlab=xlab1, ylab=ylab1, col=col1, type="l", 
         lwd=2, ...) 
  }
  
  K = length(model$lambda)
  if (is.null(col2)) {
    col2 = c("dodgerblue", "ivory4", "red4", 4:K)
  }
  x = sort(model$x)
  
  hist(x, prob=TRUE, breaks=breaks, main=main2, xlab=NA, ylab=NA, 
       axes=FALSE, ...)
  box(col="white")
  axis(side=1, tck=-0.015, mgp=c(3,0.5,0), tick=TRUE, labels=labels)
  axis(side=2, tck=-0.015, mgp=c(3,0.5,0))
  
  mtext(side=2, ylab2, line=1.5)
  if (labels == TRUE) {
    mtext(side=1, xlab2, line=1.5)
  }
  
  if(model$type == "gauss") {
    gauss = function(z) model$lambda[i] * 
      dnorm(z, mean=model$mu[i], sd=model$sigma[i]) 
    from  = x[1]
    to    = x[length(x)]

    for (i in 1:K) {
      curve(gauss, from=from, to=to, col=col2[i], lwd=2, add=TRUE)
    }
  } else if (model$type == "student"){
    student = function(z) model$lambda[i] * exp(StudentLogProb(z, 
      mu=model$mu[i], sigma=model$sigma[i], nu=model$nu[i]))
    from  = x[1]
    to    = x[length(x)]

    for (i in 1:K) {
      curve(student, from=from, to=to, col=col2[i], lwd=2, add=TRUE)
    }
  } else {
    stop("Error: unrecognized model!")
  }
}
