## Helper function for bivariate pearson correlation.
pearson = function(data) {
  cormatrix = cor(data, method="pearson")
  cormatrix[1, 2]
}

survcorr = function
(
  formula1,  # Formula for defining time-to-event 1, e.g. Surv(TIME1, STATUS1) ~ 1.
  formula2,  # Formula for defining time-to-event 2, =Surv(TIME2, STATUS2) ~ 1.
  data,    # Data set 
  methods="imi",  # Vector of method names. Allowed are "imi" and "copula", but copula not yet implemented.
  alpha=0.05,     # Confidence level; interval is [alpha, 1-alpha].
  intra=FALSE,    # TRUE if the survival times are symmetric and interchangeable.
  M=10,           # Imputation size.
  MCMCSteps=10,   # Number of MCMC steps for double-censored observations.
  epsilon= 0.001, # Precision epsilon.
  maxiter=100     # Maximum number of iterations.
)
### IMI (Iterative Multiple Imputation)
### MP, 2013-12
{
  ## Extract times ('t') and events ('delta') of data set.
  ## For <intra=T> mirror the data. 
  tmp1 = as.matrix(model.frame(formula1, data))
  if (intra) {
    tmp = as.matrix(model.frame(formula2, data))
    tmp2 = rbind(tmp, tmp1)
    tmp1 = rbind(tmp1, tmp)
  }
  else
    tmp2 = as.matrix(model.frame(formula2, data))
  
  colnames(tmp1) = colnames(tmp2) = c("time", "status")
  tmp1 = as.data.frame(tmp1)
  tmp2 = as.data.frame(tmp2)
  
  t1 = as.vector(tmp1$time)
  delta1 = as.vector(tmp1$status)
  order1 = order(-t1)
  
  t2 = as.vector(tmp2$time)
  delta2 = as.vector(tmp2$status)
  order2 = order(-t2)
  
  n = length(t1)
  minTime = 0.00000001
  
  ## Fill data set. 
  data1 = data.frame(index=1:n, t1=t1, delta1=delta1, t1unif=rep(NA, n), z1=rep(NA, n))[order1, ]
  obj = survfit(coxph(Surv(time, status) ~ 1, data=tmp1, model=TRUE), type="aalen", se.fit=FALSE, conf.type="none")
  # Get one row per object.
  data1$t1unif = rev(rep(obj$surv, obj$n.event + obj$n.censor))
  data1$t1unif = pmin(pmax(data1$t1unif, minTime), 1 - minTime)
  data1$z1 = qnorm(data1$t1unif)
  
  data2 = data.frame(index=1:n, t2=t2, delta2=delta2, t2unif=rep(NA, n), z2=rep(NA, n))[order2, ]
  obj = survfit(coxph(Surv(time, status) ~ 1, data=tmp2, model=TRUE), type="aalen", se.fit=FALSE, conf.type="none")
  # Get one row per object.
  data2$t2unif = rev(rep(obj$surv, obj$n.event + obj$n.censor))
  data2$t2unif = pmin(pmax(data2$t2unif, minTime), 1 - minTime)
  data2$z2 = qnorm(data2$t2unif)
  
  ##d = cbind(data1[order(data1$index), -1], data2[order(data2$index), -1])
  ##print(d)
  
  z1orig = data1$z1[order(data1$index)]
  z2orig = data2$z2[order(data2$index)]
  
  ## 1) Get initial correlation coefficient estimate.
  r0 = pearson(cbind(z1orig, z2orig)[delta1 & delta2, , drop=F])
  rj.all = rep(NA, M)
  
  ## Calc. indicators for censoring.
  ind1 = !delta1 & delta2
  ind2 = delta1 & !delta2
  indBoth = !delta1 & !delta2
  nBoth = sum(indBoth)
  
  z1M = z2M = matrix(NA, n, M)
  for (iImp in 1:M) {
    rj = r0
    
    ## Generate random data for 2).
    runif1 = runif(n)
    runif2 = runif(n)
    runif1mcmc = matrix(runif(nBoth*MCMCSteps), nBoth, MCMCSteps)
    runif2mcmc = matrix(runif(nBoth*MCMCSteps), nBoth, MCMCSteps)
    
    for (iter in seq(length=maxiter)) {
      ## Reset the z's
      z1 = z1orig
      z2 = z2orig
      
      ## a) Only z1i is censored.
      p1 = 1 - pnorm(z1, mean=rj * z2, sd=sqrt(1 - rj^2))     # n
      unif1 = 1 - pmin(pmax(p1 + runif1 * (1 - p1), minTime), 1 - minTime)     # n
      z1[ind1] = qnorm(unif1, mean=rj * z2, sd=sqrt(1 - rj^2))[ind1]
      
      ## b) Only z2i is censored.
      p2 = 1 - pnorm(z2, mean=rj * z1, sd=sqrt(1 - rj^2))     # n
      unif2 = 1 - pmin(pmax(p2 + runif2 * (1 - p2), minTime), 1 - minTime)     # n
      z2[ind2] = qnorm(unif2, mean=rj * z1, sd=sqrt(1 - rj^2))[ind2]
      
      ## c) Both are censored.
      z1new = z1[indBoth]
      z2new = z2[indBoth]
      for (MCMCStep in seq(length=MCMCSteps)) {
        p1 = 1 - pnorm(z1[indBoth], mean=rj * z2new, sd=sqrt(1 - rj^2))     # n
        unif1 = 1 - pmin(pmax(p1 + runif1mcmc[, MCMCStep] * (1 - p1), minTime), 1 - minTime)     # n
        z1new = qnorm(unif1, mean=rj * z2new, sd=sqrt(1 - rj^2))
        
        p2 = 1 - pnorm(z2[indBoth], mean=rj * z1new, sd=sqrt(1 - rj^2))     # n
        unif2 = 1 - pmin(pmax(p2 + runif2mcmc[, MCMCStep] * (1 - p2), minTime), 1 - minTime)     # n
        z2new = qnorm(unif2, mean=rj * z1new, sd=sqrt(1 - rj^2))
      }
      z1[indBoth] = z1new
      z2[indBoth] = z2new
      
      ## Save old correlations.
      rj.old = rj
      
      ## 3) Calculate the new correlation.
      rj = pearson(cbind(z1, z2))
      rj.all[iImp] = rj
      
      ## 5) Stop condition: all correlations converge.
      if (abs(rj - rj.old) < epsilon) {
        z1M[, iImp] = z1
        z2M[, iImp] = z2
        break
      }
    }
  }
  
  ## 6) Calculate point and interval estimates.
  n = nrow(data)  # Above n was twice as large for intra=TRUE.
  if (intra)
    df0 = n - 3/2
  else
    df0 = n - 3
  rj.trans = atanh(rj.all)
  rj.t.mean = mean(rj.trans)
  rj.t.var = 1 / df0
  BetwVar = var(rj.trans)
  TotalVar = rj.t.var + BetwVar * (M + 1) / M
  r.W.hat = tanh(rj.t.mean)
  
  df = (M - 1) * (1 + M / (BetwVar * (M + 1) * df0))^2 
  limits = tanh(rj.t.mean + c(-1, 1) * qt(1 - alpha/2, df) * sqrt(TotalVar))
  
  ## Create an object of class <survcorr>
  simData = list(M=M, z1M=z1M, z2M=z2M, delta1=delta1, delta2=delta2, t1=t1, t2=t2)
  obj = list(rho=r.W.hat, ci.lower=limits[1], ci.upper=limits[2], simData=simData, M=M, MCMCSteps=MCMCSteps, rj.trans=rj.trans, rj.t.mean=rj.t.mean, 
             var=c(within=rj.t.var, between=BetwVar, total=TotalVar), df=df, intra=intra, alpha=alpha, call=match.call())
  class(obj) = "survcorr"
  
  obj
}

print.survcorr = function
(
  x, ...
)
  ## Print the result of a <survcorr> object.
  ## MP, 2014-02
{
  print(c(rho=x$rho, lower=x$ci.lower, upper=x$ci.upper))
}

summary.survcorr = function
(
  object, ...
)
  ### Show the summary of a <survcorr> object.
  ### MP, 2014-02
{
  cat("Results of correlation analysis by Iterative Multiple Imputation\n\n")
  cat("\nContingency table of events:")
  print(table(object$simData$delta1, object$simData$delta2))
  cat("\n\nCorrelation coefficient rho and ", (1-object$alpha)*100,"% confidence interval:\n")
  print(object)
  cat("\n\nImputations: ", object$M, "\n")
  cat("MCMC steps:  ", object$MCMCSteps, "\n")
  cat("\nPosterior mean of transformed atanh(rho)\n ", object$rj.t.mean, "\n")
  cat("\nPosterior variance of atanh(rho)\n")
  print(object$var)
  cat("\nDegrees of freedom\n", object$df,"\n")

}

plot.survcorr = function
(
  x,
  what="uniform",   # Plot 'uniform' or 'times'.
  imputation=1,     # Vector of imputation indices to show. E.g. 1:10, or 1:obj$simData$M. Only used if what=='uniform'.
  xlab=switch(what, copula=expression(hat(F)(t[2])), uniform=expression(hat(F)(t[2])), times=expression(t[2])),
  ylab=switch(what, copula=expression(hat(F)(t[1])), uniform=expression(hat(F)(t[1])), times=expression(t[1])),
  xlim,             # Optional axis limits.
  ylim,             # Optional axis limits.
  main=switch(what, copula="Bivariate Copula", uniform="Bivariate Copula", times="Bivariate Survival Times"),
  legend=TRUE,      # Should appear a legend onb the upper left?
  cex.legend=switch(what, copula=0.8, uniform=0.8, times=0.7),   # Character size of plot legend.
  pch="*",          # pch value used in the plot.
  colEvent="black", # Color of events.
  colImput="gray",  # Color of imputations.
  ...
)
### Plot of a <survcorr> object (e.g. made by IMI approach).
### MP, 2014-02
{
  if(what=="copula") what<-"uniform"
  ## Extract z-values of selected imputations.
  d = x$simData
  z1M = d$z1M[, imputation, drop=FALSE]
  z2M = d$z2M[, imputation, drop=FALSE]
  nimp<-length(imputation)
  
  ## Event indicators.
  ind1 = d$delta1 & !d$delta2
  ind2 = !d$delta1 & d$delta2
  indBoth = d$delta1 & d$delta2
  indNone = !d$delta1 & !d$delta2
  
  if (what == "uniform") {
    y = 1 - pnorm(z1M)
    x = 1 - pnorm(z2M)
    if (missing(xlim)) xlim = 0:1
    if (missing(ylim)) ylim = 0:1
  }
  else if (what == "times") {
    y = d$t1
    x = d$t2
    if (missing(xlim)) xlim = range(x)
    if (missing(ylim)) ylim = range(y)
  }
  else
    stop("Value for <what> not allowed.")
  
  ## Adjust distance from axis label to axis (needed cause hat(F) was out of area)
  ## http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
  on.exit({ par(omi = c(0, 0, 0, 0)); par(mgp = c(3, 1, 0)) })
  par(oma = c(0, 1, 0, 0))
  par(mgp = c(2.2, 1, 0))  
  
  ## Draw points.
  plot(xlim, ylim, type="n", xlab=xlab, ylab=ylab, main=main, ...)
  ## Try: pch=20, cex=0.8, ...) = large dot.  ## ASCII 42.
  
  if (what == "uniform") {
    points(x[indBoth, 1], y[indBoth, 1], pch=pch, col=colEvent)
    points(x[indNone, 1:nimp], y[indNone, 1:nimp], pch=pch, col=colImput) 
    points(x[ind1], y[ind1], pch=pch, col=colImput) 
    points(x[ind2], y[ind2], pch=pch, col=colImput) 
    
    if (legend)
      legend("topleft", 
             legend=c("Event", "Censored"),
             ##pch=c(20, 158, 157, 42),  # dot, line|, line-, star
             pch=pch, # ASCII if > 32
             col=c(colEvent, colImput),
             cex=cex.legend,
             bty="n")
  }
  else if (what == "times") {
    points(x[indBoth], y[indBoth], pch=pch, col=colEvent)
    
    ## Draw small arrows. [2 alternative implementations tested]
    ##arrowLen = 0.025 * diff(ylim)
    ##arrows(x[ind1], y[ind1] - arrowLen/2, y1=y[ind1] + arrowLen/2, length=arrowLen, col=colImput)
    ##arrowLen = 0.025 * diff(xlim)
    ##arrows(x[ind2] - arrowLen/2, y[ind2], x1=x[ind2] + arrowLen/2, length=arrowLen, col=colImput)
    arrowLen = 0.02
    length = 0.05
    arrow.plot(x[ind1], y[ind1], u=1, v=0, arrow.ex=arrowLen, length=length, col=colImput) 
    arrow.plot(x[ind2], y[ind2], u=0, v=1, arrow.ex=arrowLen, length=length, col=colImput)
    arrow.plot(x[indNone], y[indNone], u=1, v=1, arrow.ex=arrowLen, length=length, col=colImput)
    
    if (legend)
      legend("topright", 
             legend=c("Uncensored times 1&2", 
                      "Censored time 1", 
                      "Censored time 2", 
                      "Censored times 1&2"), 
             pch=c(pch, "|", "-", "/"),
             col=c(colEvent, colImput, colImput, colImput),
             cex=cex.legend,
             bty="n")
  }
}


