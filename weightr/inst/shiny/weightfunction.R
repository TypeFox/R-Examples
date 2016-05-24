weightfunction <- function(effect, v, npred, steps, XX, prednames, weights=NULL) {

  si <- sqrt(v)
  effect <- effect
  number <- length(effect)
  v <- v
  p <- 1-pnorm(effect/sqrt(v))

  neglike1 <- function(pars) {
    vc = pars[1]
    beta = pars[2:(npred+2)]
    mn = XX%*%beta
    eta = sqrt(v + vc)
    b = 1/2 * sum(((effect-mn)/eta)^2)
    c = sum(log(eta))
    return(b+c)
  }

  neglike2 <- function(pars) {
    vc = pars[1]
    beta = pars[2:(npred+2)]
    if(is.null(weights)==FALSE){
      w <- weights
    }
    else{
    w = c(1,pars[(npred+3):length(pars)])
    }
    mn = XX%*%beta
    a = sum(log(w[wt]))
    eta = sqrt(v + vc)
    b = 1/2 * sum(((effect-mn)/eta)^2)
    c = sum(log(eta))
    Bij <- matrix(rep(0,number*nsteps),nrow=number,ncol=nsteps)
    bi = -si * qnorm(steps[1])
    Bij[,1] = 1-pnorm((bi-mn)/eta)
    if(nsteps > 2){
      for(j in 2:(length(steps)-1)) {
        bi = -si * qnorm(steps[j])
        bilast = -si * qnorm(steps[j-1])
        Bij[,j] = pnorm((bilast-mn)/eta) - pnorm((bi-mn)/eta)
      }
    }
    bilast = -si * qnorm(steps[length(steps)-1])
    Bij[,length(steps)] = pnorm((bilast-mn)/eta)

    swbij = 0
    for(j in 1:length(steps)) swbij = swbij + w[j]*Bij[,j]
    d = sum(log(swbij))
    # (Uncomment if needed to see what's going wrong.)  print(pars)
    return(-a + b + c + d)
  }

  gradient1 <- function(pars) {
    vc = pars[1]
    beta = pars[2:(npred+2)]
    mn = XX%*%beta
    eta2 = v + vc
    a = 0.5*sum(1/eta2 - (effect-mn)^2/eta2^2)
    b = rep(0,(1+npred))
    for(i in 1:length(b)) b[i] = sum( (effect-mn)*-XX[,i]/eta2 )
    return(matrix(c(a,b),nrow=1+length(b),ncol=1))
  }

  intervaltally <- function(p, steps) {
    p1 <- cut(p, breaks=c(-Inf,steps), labels=steps)
    return(p1)
  }

  #UNADJUSTED WITHOUT MODERATORS
  if(steps == 600 && npred == 0 && XX == 600){

    nsteps <- 0
    XX <- cbind(rep(1,number))
    pars <- c(mean(v)/4,mean(effect))

    output1 <- optim(par=pars,fn=neglike1, lower=c(0,rep(-Inf,(npred+1))),method="L-BFGS-B",hessian=TRUE)

    if(output1$convergence == 1){
      print("Maximum iterations reached without convergence. Consider re-examining your model.")
    }
    if(output1$convergence >= 51){
      print("The model has failed to converge and produced an error. Consider re-examining your model.")
    }

    Parameters <- output1$par
    Standard_Errors <- sqrt(diag(solve(output1$hessian)))

    results <- cbind(c(output1$par[1], output1$par[2]),Standard_Errors)

    #results <- cbind(c("Parameters",output1$par[1], output1$par[2]),c("Standard Errors", sqrt(diag(solve(output1$hessian)))))
    #return(grid.table(results, gp=gpar(fontsize=20), rows=c("Variance Component", "Intercept"), cols=c("Estimate", "Standard Error")))
    resultsb <- data.frame(results, row.names=c("Variance Component", "Intercept"))
    colnames(resultsb) <- c("Parameters", "Standard Errors")
    return(resultsb)

  }

  #ADJUSTED WITHOUT MODERATORS
  if(steps != 600 && npred == 0 && XX == 600){



    XX <- cbind(rep(1,number))
    nsteps <- length(steps)
    pars <- c(mean(v)/4, mean(effect), rep(1,nsteps-1))

    wt <- rep(1,number)
    for(i in 1:number) {
      for(j in 2:nsteps) {
        if (-si[i]*qnorm(steps[j]) <= effect[i] && effect[i] <= -si[i]*qnorm(steps[j-1])) wt[i] = j
      }
      if(  effect[i] <= -si[i]*qnorm(steps[nsteps-1])) wt[i] = nsteps
    }

    output2 <- optim(par=pars,fn=neglike2,
                     lower=c(0,rep(-Inf,1),rep(0.01,(nsteps-1))),
                     method="L-BFGS-B",hessian=TRUE)

    if(output2$convergence == 1){
      print("Maximum iterations reached without convergence. Consider re-examining your model.")
    }
    if(output2$convergence >= 51){
      print("The model has failed to converge and produced an error. Consider re-examining your model.")
    }

    Parameters <- output2$par

    results <- matrix(nrow=(length(output2$par)),ncol=2)

    if(is.null(weights)){
    results[,1] <- output2$par
    results[,2] <- sqrt(diag(solve(output2$hessian)))
    }
    else{
      results[,1] <- c(output2$par[1:2],weights[2:length(weights)])
      results[,2] <- rep(NA, length(output2$par))
    }
    rowlabels <- c(0, length(output2$par))
    rowlabels[1] <- "Variance Component"
    rowlabels[2] <- "Intercept"
    for(i in 3:(nsteps + 1)){
      rowlabels[i] <- paste(c(steps[i - 2], "< p-values <", steps[i - 1], "weight"), collapse=" ")
    }
    resultsb <- data.frame(results, row.names=c(rowlabels))
    colnames(resultsb) <- c("Parameters", "Standard Errors")
    #return(grid.table(results, gp=gpar(fontsize=20), rows=c(rowlabels), cols=c("Estimate", "Standard Error")))
    return(resultsb)
  }

  #UNADJUSTED WITH MODERATORS
  if(steps == 600 && npred >= 0 && XX != 600){

    nsteps <- 0
    pars <- c(mean(v)/4,mean(effect),rep(0,npred))

    output1 <- optim(par=pars,fn=neglike1, lower=c(0,rep(-Inf,(npred+1))),method="L-BFGS-B",hessian=TRUE)

    if(output1$convergence == 1){
      print("Maximum iterations reached without convergence. Consider re-examining your model.")
    }
    if(output1$convergence >= 51){
      print("The model has failed to converge and produced an error. Consider re-examining your model.")
    }

    results <- matrix(nrow=(length(output1$par)),ncol=2)

    results[,1] <- output1$par
    results[,2] <- sqrt(diag(solve(output1$hessian)))
    predictornumber <- seq(1,npred,by=1)
    rowlabels <- c(0, length(output1$par))
    rowlabels[1] <- "Variance Component"
    rowlabels[2] <- "Intercept"
    for(i in 3:(npred + 2)){
      rowlabels[i] <- prednames[i - 1]
    }
    resultsb <- data.frame(results, row.names=c(rowlabels))
    colnames(resultsb) <- c("Parameters", "Standard Errors")
    #return(grid.table(results, gp=gpar(fontsize=20), rows=c(rowlabels), cols=c("Estimate", "Standard Error")))
    return(resultsb)
  }

  #ADJUSTED WITH MODERATORS

  if(steps != 600 && npred >= 0 && XX != 600){

    nsteps <- length(steps)
    pars <- c(mean(v)/4,mean(effect),rep(0,npred), rep(1,nsteps-1))
    wt <- rep(1,number)
    for(i in 1:number) {
      for(j in 2:nsteps) {
        if (-si[i]*qnorm(steps[j]) <= effect[i] && effect[i] <= -si[i]*qnorm(steps[j-1])) wt[i] = j
      }
      if(  effect[i] <= -si[i]*qnorm(steps[nsteps-1])) wt[i] = nsteps
    }
    output2 <- optim(par=pars,fn=neglike2,
                     lower=c(0,rep(-Inf,(npred+1)),rep(0.01,(nsteps-1))),
                     method="L-BFGS-B",hessian=TRUE)

    if(output2$convergence == 1){
      print("Maximum iterations reached without convergence. Consider re-examining your model.")
    }
    if(output2$convergence >= 51){
      print("The model has failed to converge and produced an error. Consider re-examining your model.")
    }

    results <- matrix(nrow=(length(output2$par)),ncol=2)

    if(is.null(weights)){
      results[,1] <- output2$par
      results[,2] <- sqrt(diag(solve(output2$hessian)))
    }
    else{
      results[,1] <- c(output2$par[1:(npred+2)],weights[2:length(weights)])
      results[,2] <- rep(NA, length(output2$par))
    }
    predictornumber <- seq(1,npred,by=1)
    rowlabels <- c(0, length(output2$par))
    rowlabels[1] <- "Variance Component"
    rowlabels[2] <- "Intercept"
    for(i in 3:(npred + 2)){
      rowlabels[i] <- prednames[i - 1]
    }
#     for(i in 3:(npred + 2)){
#       rowlabels[i] <- paste(c("Slope #", predictornumber[i - 2]), collapse=" ")
#     }
    for(i in (3+npred):(length(output2$par))){
      rowlabels[i] <- paste(c(steps[i - (2 + npred)], "< p-values <", steps[i - (1 + npred)], "weight"), collapse=" ")
    }
    resultsb <- data.frame(results, row.names=c(rowlabels))
    colnames(resultsb) <- c("Parameters", "Standard Errors")
   # return(grid.table(results, gp=gpar(fontsize=20), rows=c(rowlabels), cols=c("Estimate", "Standard Error")))
   return(resultsb)

  }


  }

likelihoodfunct <- function(effect, v, npred, steps, XX) {

  si <- sqrt(v)
  effect <- effect
  number <- length(effect)
  v <- v
  p <- 1-pnorm(effect/sqrt(v))

  neglike1 <- function(pars) {
    vc = pars[1]
    beta = pars[2:(npred+2)]
    mn = XX%*%beta
    eta = sqrt(v + vc)
    b = 1/2 * sum(((effect-mn)/eta)^2)
    c = sum(log(eta))
    return(b+c)
  }

  neglike2 <- function(pars) {
    vc = pars[1]
    beta = pars[2:(npred+2)]
    w = c(1,pars[(npred+3):length(pars)])
    mn = XX%*%beta
    a = sum(log(w[wt]))
    eta = sqrt(v + vc)
    b = 1/2 * sum(((effect-mn)/eta)^2)
    c = sum(log(eta))
    Bij <- matrix(rep(0,number*nsteps),nrow=number,ncol=nsteps)
    bi = -si * qnorm(steps[1])
    Bij[,1] = 1-pnorm((bi-mn)/eta)
    if(nsteps > 2){
      for(j in 2:(length(steps)-1)) {
        bi = -si * qnorm(steps[j])
        bilast = -si * qnorm(steps[j-1])
        Bij[,j] = pnorm((bilast-mn)/eta) - pnorm((bi-mn)/eta)
      }
    }
    bilast = -si * qnorm(steps[length(steps)-1])
    Bij[,length(steps)] = pnorm((bilast-mn)/eta)

    swbij = 0
    for(j in 1:length(steps)) swbij = swbij + w[j]*Bij[,j]
    d = sum(log(swbij))
    # (Uncomment if needed to see what's going wrong.)  print(pars)
    return(-a + b + c + d)
  }

  gradient1 <- function(pars) {
    vc = pars[1]
    beta = pars[2:(npred+2)]
    mn = XX%*%beta
    eta2 = v + vc
    a = 0.5*sum(1/eta2 - (effect-mn)^2/eta2^2)
    b = rep(0,(1+npred))
    for(i in 1:length(b)) b[i] = sum( (effect-mn)*-XX[,i]/eta2 )
    return(matrix(c(a,b),nrow=1+length(b),ncol=1))
  }

  intervaltally <- function(p, steps) {
    p1 <- cut(p, breaks=c(-Inf,steps), labels=steps)
    return(p1)
  }

  #WITHOUT MODERATORS
  if(npred == 0 && XX == 600){

    nsteps <- 0
    XX <- cbind(rep(1,number))
    pars <- c(mean(v)/4,mean(effect))

    output1 <- optim(par=pars,fn=neglike1, gr=gradient1, lower=c(0,rep(-Inf,(npred+1))),method="L-BFGS-B",hessian=TRUE)

    nsteps <- length(steps)
    pars <- c(mean(v)/4, mean(effect), rep(1,nsteps-1))

    wt <- rep(1,number)
    for(i in 1:number) {
      for(j in 2:nsteps) {
        if (-si[i]*qnorm(steps[j]) <= effect[i] && effect[i] <= -si[i]*qnorm(steps[j-1])) wt[i] = j
      }
      if(  effect[i] <= -si[i]*qnorm(steps[nsteps-1])) wt[i] = nsteps
    }

    output2 <- optim(par=pars,fn=neglike2,
                     lower=c(0,rep(-Inf,1),rep(0.01,(nsteps-1))),
                     method="L-BFGS-B",hessian=TRUE)
    lrchisq <- 2*(output1$value - output2$value)
    results <- cbind(output1$value, output2$value, lrchisq, (nsteps-1), (1-pchisq(lrchisq,(nsteps-1))) )
    likelihood <- data.frame(results, row.names=c("Estimate"))
    colnames(likelihood) <- c("Unadjusted Likelihood","Adjusted Likelihood", "2*Difference","df", "p-value")

    return(likelihood)
  }

  #UNADJUSTED WITH MODERATORS
  if(npred >= 0 && XX != 600){

    nsteps <- 0
    pars <- c(mean(v)/4,mean(effect),rep(0,npred))

    output1 <- optim(par=pars,fn=neglike1, gr=gradient1, lower=c(0,rep(-Inf,(npred+1))),method="L-BFGS-B",hessian=TRUE)

    nsteps <- length(steps)
    pars <- c(mean(v)/4,mean(effect),rep(0,npred), rep(1,nsteps-1))
    wt <- rep(1,number)
    for(i in 1:number) {
      for(j in 2:nsteps) {
        if (-si[i]*qnorm(steps[j]) <= effect[i] && effect[i] <= -si[i]*qnorm(steps[j-1])) wt[i] = j
      }
      if(  effect[i] <= -si[i]*qnorm(steps[nsteps-1])) wt[i] = nsteps
    }
    output2 <- optim(par=pars,fn=neglike2,
                     lower=c(0,rep(-Inf,(npred+1)),rep(0.01,(nsteps-1))),
                     method="L-BFGS-B",hessian=TRUE)

    lrchisq <- 2*(output1$value - output2$value)
    results <- cbind(output1$value, output2$value, lrchisq, (nsteps-1), (1-pchisq(lrchisq,(nsteps-1))) )
    likelihood <- data.frame(results, row.names=c("Estimate"))
    colnames(likelihood) <- c("Unadjusted Likelihood","Adjusted Likelihood", "2*Difference","df", "p-value")
    return(likelihood)

  }

}


sampletable <- function(p, pvalues, steps){
  nsteps <- length(steps)
  results <- matrix(nrow=length(pvalues),ncol=1)
  results[,1] <- pvalues
  rowlabels <- c(0, length(results[,1]))
  rowlabels[1] <- paste(c("p-values <", steps[1]), collapse="")
  for(i in 2:nsteps){
    rowlabels[i] <- paste(c(steps[i - 1], "< p-values <", steps[i]), collapse=" ")
  }
  resultsb <- data.frame(results, row.names=c(rowlabels))
  colnames(resultsb) <- c("Number of Effects")
  return(resultsb)
}
