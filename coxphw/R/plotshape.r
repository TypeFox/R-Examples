plotshape <- function(obj, plot.x = NA, ref = NA, plot = TRUE, variable = NULL, treatment = NULL, ref.type = "value", 
                      addci = TRUE, exp = FALSE, xlab = NULL, ylab = NULL, 
                      pval = FALSE, digits = 4, ...)
{
  
  ### version 2014-02-25:     
  ### version 2010-04-29: NOT COMPATIBLE TO OLDER VERSIONS!!!
    
    ### obj:   coxphf.beta object
    ### plot.x: range of x axis
    ### ref:    reference value (only for ref.type=="value" mode)
    ### plot:   if function should be plotted
    ### treatment: variable which is in interaction with fp variable (use "")
    ### ref.type:  "value" for simple fp term (no interaction),
    ###            "interaction.time" for interaction for treatment with fp(time)
    ###            "interaction.treat" for interaction of treatment with fp(variable)
    ###            any value to specify the level of the treatment variable for which the fp of "variable" should be plotted
    ### variable:  name of fp variable (use "")

    data<-obj$dataline      # a dataline from original data set of analysis
    if (is.null(variable)) stop("No variable specified.\n")
    if (ref.type=="interaction.time" & is.null(treatment)) stop("No interaction variable for time specified.\n")
    if (is.na(ref)) ref <- min(plot.x)
    if (any(plot.x <= 0) & ref.type=="interaction.time") warning(paste("model.frame() might cancel lines with ", variable, "<=0..."))
    n.x <- length(plot.x)
    data <- data[1,]
    data[,variable]<-NA
    if(ref.type=="interaction.time" | ref.type=="interaction.treat") data[,treatment]<-1
    else if(ref.type != "value") data[,treatment]<-ref.type
    i.x <- which(is.na(data))
    plotData <- data[rep(1, n.x), , drop = FALSE]
    plotData[, i.x] <- plot.x
    rownames(plotData) <- 1:n.x
    objP <- decomposeSurv(obj$formula, plotData, sort = FALSE, obj$offset, PTpreset = obj$PTcoefs)
    kk <- ncol(objP$mm1)
    objP$mm1 <- objP$mm1[, obj$ind[1:kk], drop = FALSE]
#    objP$covnames <- objP$covnames[obj$ind]
    objP$timedata <- objP$timedata[, obj$ind[-(1:kk)], drop = FALSE]
#    plot.y <- (cbind(objP$mm1, objP$timedata)) %*% obj$coefficients
    
    x1 <- cbind(objP$mm1, objP$timedata)
   
    if (ref.type != "interaction.time" & ref.type != "interaction.treat"){
        plotDataRef <- plotData
        plotDataRef[, i.x] <- ref
        objPRef <- decomposeSurv(obj$formula, plotDataRef, sort = FALSE,
        obj$offset, PTpreset = obj$PTcoefs)
        objPRef$mm1 <- objPRef$mm1[, obj$ind[1:kk], drop = FALSE]
#        objPRef$covnames <- objPRef$covnames[obj$ind]
        objPRef$timedata <- objPRef$timedata[, obj$ind[-(1:kk)], drop = FALSE]
#        plot.yRef <- (cbind(objPRef$mm1, objPRef$timedata)) %*% obj$coefficients
#        plot.y <- plot.y - plot.yRef
        
        x0 <- cbind(objPRef$mm1, objPRef$timedata)
        diff <- x1-x0
        }
   if (ref.type == "interaction.treat"){
        plotDataRef <- plotData
        plotDataRef[, treatment] <- 0
        objPRef <- decomposeSurv(obj$formula, plotDataRef, sort = FALSE,
        obj$offset, PTpreset = obj$PTcoefs)
        objPRef$mm1 <- objPRef$mm1[, obj$ind[1:kk], drop = FALSE]
#        objPRef$covnames <- objPRef$covnames[obj$ind]
        objPRef$timedata <- objPRef$timedata[, obj$ind[-(1:kk)], drop = FALSE]
#        plot.yRef <- (cbind(objPRef$mm1, objPRef$timedata)) %*% obj$coefficients
#        plot.y <- plot.y - plot.yRef
        
        x0 <- cbind(objPRef$mm1, objPRef$timedata)
        diff <- x1-x0
   }

   if (ref.type == "interaction.time") { diff <- x1 }    

   plot.y <- diff %*% obj$coefficients
   se2 <- sapply(1:nrow(diff), function(I) t(diff[I,]) %*% obj$var %*% diff[I,])                
   gammavt <- sqrt(qchisq(1-obj$alpha, df=1) * se2)                                             
   cilower <- plot.y - gammavt
   ciupper <- plot.y + gammavt

   if (pval & !plot) { p <- 1 - pchisq(plot.y^2 / se2, df = 1) } else { p <- NA }
    
   refline <- 0
   if (exp) {
     plot.y <- exp(plot.y)
     cilower <- exp(cilower)
     ciupper <- exp(ciupper)
     if(is.null(ylab)) { ylab <- "relative hazard" }
     refline <- 1
   }
        
   if (plot) {
        if(is.null(xlab)) { xlab <- variable }
        if(is.null(ylab)) { ylab <- "log relative hazard" }
        if (addci) {ylimit <- range(plot.y, ciupper, cilower) } else {ylimit <- range(plot.y) }
        plot(plot.x, plot.y, lty = 1, xlab=xlab, ylab=ylab, ylim=ylimit, ...)
        if (addci) {
          lines(x=plot.x, y=cilower, lty=2)
          lines(x=plot.x, y=ciupper, lty=2) 
        }
        abline(refline, 0, col="gray", lty = 2)
    }
    
    
    object <- list(xbeta=plot.y, ci.lower=cilower, ci.upper=ciupper, p=p, alpha=obj$alpha, plot.x=plot.x, exp=exp)    
    
    if (!plot)  {
      printobj <- data.frame(plot.x, plot.y, cilower, ciupper)
      if (exp)  { dimnames(printobj)[[2]] <- c(variable, "HR", paste(c("HR lower", "HR upper"), 1-obj$alpha, sep=" ")) } else
      if (!exp) { dimnames(printobj)[[2]] <- c(variable, "coef", paste(c("coef lower", "coef upper"), 1-obj$alpha, sep=" ")) }
      
      if (pval) { printobj <- cbind(printobj, p) }
      print(round(printobj, digits = digits))
    }
    
    return(invisible(object))
}
