`relaxo` <-
function(X, Y, phi= seq(0,1, length=4), max.steps= min(2*length(Y),2*ncol(X)),
                   fast = TRUE,  keep.data = TRUE, warn=TRUE )
{
  Y <- as.numeric(Y)
  if(warn){
    if( abs(mean(Y))> 0.01*sd(Y)) warning("response variable not centered")
    if( any( abs(apply(X,2,mean)) > 0.01* apply(X,2,sd) )) warning("predictor variables not centered")
    if( sd(as.numeric(apply(X,2,sd)))>0.001) warning("predictor variables not scaled")
  }
  
  stopifnot(require("lars")) # bail out if package "lars" cannot be loaded

  if(max(phi) > 1 || min(phi) < 0) stop("phi must be in [0,1]")
  if(!any(phi == 1)){
    warning("adding 1 to phi (Lasso solution)")
    phi <- c(phi,1)
  }
  phi <- sort(phi, decreasing=TRUE) # ==>  phi[1] == 1

  stopifnot((n.phi <- length(phi)) >= 2) 



  larsfit <- lars(X, Y, use.Gram=FALSE, max.steps=max.steps)
  
  beta <- coef(larsfit)
  
  nSteps <- nrow(beta)
  if(nSteps < 3)
      stop("'Need at least 3 lars() steps; maybe increase 'max.steps'?")
  betarelaxo <- matrix(nrow = n.phi*(nSteps-1), ncol=ncol(beta))

  lambda <- numeric(nSteps)
  for(s in 2:nSteps)
    lambda[s] <- abs(crossprod(Y - X %*% beta[s,], X[, which(beta[s,] != 0)[1]]))
  lambda[1] <- max(abs( crossprod(Y, X )))


  lambdaall <- rep(lambda[1:(length(lambda)-1)], each=n.phi)
  phiall <- rep(phi, nSteps-1)

  betarelaxo[1:n.phi, ] <- rep(0, ncol(beta))

  for (s in 2:(nSteps-1)) {

    betarelaxo[ n.phi*(s-1) + 1, ] <- beta[s, ]

    diffbeta <- beta[s, ] - beta[s-1, ]
    difflambda <- lambda[s-1] - lambda[s]

    ## try shooting method
    OLStrybeta <- beta[s, ] + diffbeta* lambda[s]/difflambda

    if( !all(sign(OLStrybeta) == sign(beta[s, ])) && !fast ) {
      ## new estimation necessary
      Xreduced <- X[, beta[s, ] != 0]
      larsnew <- lars(Xreduced, Y, use.Gram=FALSE)
      betanew <- coef(larsnew)
      nStepsnew <- nrow(betanew)
      lambdanew <- rep(0, nStepsnew)
      lambdanew[1] <- max( abs( crossprod(Y, Xreduced) ) )
      if(nStepsnew >= 2) for (snew in 2:nStepsnew) {
        ## take the first of the selected variables
        selVar <-  which(betanew[snew, ] != 0)[1]
        lambdanew[snew] <- abs( crossprod(Y - Xreduced %*% betanew[snew, ],
                                          Xreduced[, selVar] ))
      }
      phinew <- lambdanew/max(lambdanew)

      for (phic in 2:n.phi) {

        select <- which.min( (phi[phic] - phinew)^2 )

        relaxed <- rep(0, ncol(X))
        relaxed[ beta[s, ] != 0 ] <- coef(larsnew, s=select, mode="step")
        betarelaxo[n.phi*(s-1) + phic , ] <- relaxed
      }

    } else {
      ## shooting method sucessful
      for (phic in 2:n.phi) {
        betarelaxo[n.phi*(s-1) + phic , ] <-
            phi[phic]*beta[s, ] + (1-phi[phic])* OLStrybeta
      }
    }

  }

  relaxo <- c(if(keep.data) list(X = X, Y = Y), 
              list(beta = betarelaxo,
                   lambda = lambdaall,
                   phi = phiall))
  class(relaxo) <- "relaxo"

  return(relaxo)
}

