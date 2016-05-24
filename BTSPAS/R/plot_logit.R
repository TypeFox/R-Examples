# 2014-09-01  CJS First edition of this function
# Take the input values and create a ggplot object for the logitP's with the credible intervals plotted
# Input are the usual data values along with the MCMC results



plot_logitP <- function(title, time, n1, m2, u2, logitP.cov, results){
  #  Plot the observed and fitted logit(p) values along with posterior limits 
  #  n1, m2, u2 are the raw data (u2 has been adjusted upward for sampling fraction < 1 prior to call)
  #  logitP.cov is the covariate matrix for modelling the logit(P)'s
  #  results is the summary table from WinBugs
  #

  Nstrata.rel <- length(n1)
  Nstrata.cap <- length(u2) 

  # which rows of the result summary contain the logitP[xx] ?
  results.row.names <- rownames(results$summary)
  logitP.row.index  <- grep("^logitP", results.row.names)
  logitP.res<- as.data.frame(results$summary[logitP.row.index,])  # summary statistics 
  
  # We need to extract the time index from the row.names
  logitP.res$time.index <- aaply(rownames(logitP.res), 1, function(x){
     # extract the time index from logitP[xx]
     temp <- unlist(strsplit(x,    "[", fixed=TRUE))[2]
     temp <- unlist(strsplit(temp, "]", fixed=TRUE))[1]
     temp <- as.numeric(temp)
     temp
  })
  
  # Only retain entries in the range of 1... length(u2)
  logitP.res <- logitP.res[logitP.res$time.index <= Nstrata.cap,]
  logitP.res$time <- time[logitP.res$time.index]
  
  # Set up the bottom axis title
  xtitle <- paste("Time\nHorizontal line is estimated beta.logitP[1]",
                     "\nInner fence is c.i. on beta.logitP[1]",
                     "\nOuter fence is 95% range on logit(p)")
  if(ncol(as.matrix(logitP.cov))>1){
     xtitle<-paste(xtitle,"\nDashed line is second covariate")}

  # Extract the upper and lower ci
  logitP.res$lcl <- logitP.res[, "2.5%"]
  logitP.res$ucl <- logitP.res[,"97.5%"]
  
  myplot <- ggplot(data=logitP.res, aes(x=time, y=mean))+
    ggtitle( paste(title,"\nPlot of logit(p[i]) with 95% credible intervals"))+
    xlab(xtitle)+ylab("logit(p) + 95% credible interval")+
    geom_point(size=3)+
    geom_line()+
    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)

  # If this is a non-diagonal case, also plot the raw logits
  if(!is.matrix(m2)){
    raw_logitP <- logit((m2+1)/(n1+2))
    myplot <- myplot + annotate("point", x=time, y=raw_logitP, shape=1) 
  }        # based on raw data
  
  # plot the posterior mean of the logitP if there is only one column for a covariate   
  if(ncol(as.matrix(logitP.cov))==1){  # if only 1 column for covariate vector, usually an intercept
    # plot the posterior mean of the beta.logitP[1] term which is usually
    #      the intercept in most models with covariates along with 95% credible interval
    intercept.row.index    <- grep("beta.logitP[1]", results.row.names, fixed=TRUE)
    intercept <- results$summary[intercept.row.index,]
    mean<- intercept["mean"]
    lcl <- intercept["2.5%"]
    ucl <- intercept["97.5%"]
    
    myplot <- myplot +
        geom_hline(yintercept=mean)+
        geom_hline(yintercept=lcl, linetype=2)+
        geom_hline(yintercept=ucl, linetype=2)
     
    # plot the posterior "95% range" for the logit(P)'s based on N(xip, sigmaP^2)
    sigmaP.row.index <- grep("sigmaP", results.row.names)
    sigmaP <- results$summary[sigmaP.row.index,]
    lcl <- intercept["mean"]-2*sigmaP["mean"]
    ucl <- intercept["mean"]+2*sigmaP["mean"]
    myplot <- myplot +
      geom_hline(yintercept=lcl, linetype=3)+
      geom_hline(yintercept=ucl, linetype=3)
  }
  
  # plot residuals of the logit(P)'s against the various covariates
  # to be done in my next life
  return(myplot)
}
