##IterGS.r
###This function conducts an iterative graphlet screening.


IterGS<-function(y.tilde, gram, gram.bias, cg.all, sp, tau, nm, q0=0.1, scale=1, max.iter=3, std.thresh=1.05, beta.initial=NULL) {
  ########################################################################
  # Run the Graphlet Screening procedure and do the refinement up to max.iter times. 
  #
  # Args:
  #   y.tilde: X'Y
  #   gram: the threholded gram matrix
  #   gram.bias: the bias of the threholded gram matrix
  #   cg.all: all the connected cg.alls of gram with size no more than nm.
  #   sp: the expected sparse level
  #   tau: the minimal signal strength to be detected
  #   nm: the maximal size of the connected subgaphs considered in the screening step.
  #   q0: the minimal screening parameter.
  #   scale: optional numerical parameter of the screening step. The default is 1
  #   max.iter: the maximal number of iterations. The default is 3.
  #   std.thresh: the threshold of the std change that stop the loop. The default is 1.05.
  #   beta.initial: the initial estimate of beta in reducing the bias. 
  #               The default is uu*sign(y.tilde)*(abs(y.tilde)>uu).
  #  
  # Returns:
  #   estimate: the GS estimate for beta
  #   n.iter: the number o iterations it uses
  #
  ###############################################################################

  p <- dim(gram)[2]
  r <- tau^2/2/log(p)
  v <- 1-log(sp)/log(p)
  uu <- sqrt(2*r*log(p))
  lambda <- sqrt(2*v*log(p))

  # use a conservative hard thresholding result as the initial estimate
  # such arrangement helps reduce the bias.
  if(length(beta.initial) == 0) {
    beta.initial <- uu*sign(y.tilde)*(abs(y.tilde)>uu)
  }
  beta.gs <- beta.initial

  # refinement interations
  w <- y.tilde # initial

  for (it in 1:max.iter){          
    # prepare  
    last.w.std <- sqrt((sum(w^2)-p*mean(w)^2)/(p-1))
    last.beta <- beta.gs
    last.nonzero <- (last.beta != 0)
    # update w
    w <- y.tilde-gram.bias[, last.nonzero] %*% last.beta[last.nonzero]
    if (sqrt((sum(w^2)-p*mean(w)^2)/(p-1)) > std.thresh * last.w.std) { 
      ## abort the current iteratoin
      n.iteration <- it - 1
      break 
    }
    ## update beta.gs
    survivor <- ScreeningStep(w, gram, cg.all, nm, v, r, q0, scale)
    beta.gs <- CleaningStep(survivor, w, gram, lambda, uu)
    n.iteration <- it ## count finished iterations         
  } 
  return(list(estimate=beta.gs, n.iter=n.iteration))
}

