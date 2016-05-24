balance <- function (rawdata, matched, pscore.fit, factor=TRUE)
{

    #int <- attr(terms(pscore.fit), "intercept")
    call.raw <- call.matched <- pscore.fit$call
    call.raw$data <- substitute(rawdata)
    call.matched$data <- substitute(matched)
    if(!is.call(pscore.fit$call$formula)){
        call.raw$formula <- formula(terms(pscore.fit))
    }

    if (!factor){
        form <- gsub("factor(", "", call.raw$formula, fixed = TRUE)
        form <- gsub(")", "", form, fixed = TRUE)
        form <- as.formula(paste(form[2], form[1], form[3], -1))
        call.raw$formula <- call.matched$formula <- form
    }
    else{
        form <- call.raw$formula
        form <- as.formula(paste(form[2], form[1], form[3], -1))
        call.raw$formula <- call.matched$formula <- form
    }
    fit.raw <- eval(call.raw)
    fit.matched <- eval(call.matched)
    class(fit.raw$formula) <- class(fit.matched$formula) <- c("bayesglm", "formula")
    treat.raw <- fit.raw$y
    treat.matched <- fit.matched$y
    pred.raw <- model.matrixBayes(fit.raw$formula, data=rawdata, keep.order=TRUE)
    pred.matched <- model.matrixBayes(fit.matched$formula, data=matched, keep.order=TRUE)

    #if (int){
    #    pred.raw <- model.matrix(fit.raw)[,-1]
    #    pred.matched <- model.matrix(fit.matched)[,-1]
    #}
    #if (!int){
    #    pred.raw <- model.matrix(fit.raw)
    #    pred.matched <- model.matrix(fit.matched)
    #}

    if(dim(pred.raw)[2]!=dim(pred.matched)[2])
        warnings("number of covariates of the raw data does not equal to
            that of the matched data! This might be due to the drop of
            factor levels.  Use factor=FALSE to proceed!")

    raw.dat <- data.frame(pred.raw,treat=treat.raw)
    matched.dat <- data.frame(pred.matched,treat=treat.matched)
    covnames <- c(colnames(pred.matched),"treat")
    K <- length(covnames)-1

    # diff.mean.rawdata
    diff.means=matrix(NA,K,6)
    for(i in 1:K){
        diff.means[i,1:2] <- c(mean(raw.dat[(raw.dat[,"treat"]==1),i]),
            mean(raw.dat[(raw.dat[,"treat"]==0),i]))
        diff.means[i,3] <- diff.means[i,1]-diff.means[i,2]
        diff.means[i,5] <- sqrt(var(raw.dat[(raw.dat[,"treat"]==1),i])/
            sum((raw.dat[,"treat"]==1)) + var(raw.dat[(raw.dat[,"treat"]==0),i])/sum((raw.dat[,"treat"]==0)))
        diff.means[i,6] <- sqrt((var(raw.dat[(raw.dat[,"treat"]==1),i])+
            var(raw.dat[(raw.dat[,"treat"]==0),i]))/2)
        diff.means[i,4] <- diff.means[i,3]/diff.means[i,6]
    }
    dimnames(diff.means) <- list(covnames[-(K+1)],
        c("Treat","control","diff","diff.std","se","sd"))

    # diff.means.matched.dat
    diff.means.matched=matrix(NA,K,6)
    for(i in 1:K){
        diff.means.matched[i,1:2] <- c(mean(matched.dat[(matched.dat[,"treat"]==1),i]),
            mean(matched.dat[(matched.dat[,"treat"]==0),i]))
        diff.means.matched[i,3] <- diff.means.matched[i,1]-diff.means.matched[i,2]
        diff.means.matched[i,5] <- sqrt(var(matched.dat[(matched.dat[,"treat"]==1),i])/
            sum((raw.dat[,"treat"]==1)) + var(raw.dat[(raw.dat[,"treat"]==0),i])/sum((raw.dat[,"treat"]==0)))
        diff.means.matched[i,6] <- sqrt((var(raw.dat[(raw.dat[,"treat"]==1),i])+
            var(raw.dat[(raw.dat[,"treat"]==0),i]))/2)
        diff.means.matched[i,4] <- diff.means.matched[i,3]/diff.means.matched[i,6]
    }
    dimnames(diff.means.matched) <- list(covnames[-(K+1)],
        c("Treat","control","diff","diff.std","se","sd"))

    out <- list(diff.means.raw=diff.means,
      diff.means.matched=diff.means.matched, covnames=covnames)
    class(out) <- "balance"
    return(out)
}


print.balance <- function(x, ..., digits= 2)
{
 cat("Differences in Means of Unmatched Data\n")
 cat("--\n")
 print(round(x$diff.means.raw, digits=digits))
 cat("--\n")
 cat("\n")
 cat("Differences in Means of Matched Data\n")
 cat("--\n")
 print(round(x$diff.means.matched, digits=digits))
 cat("--\n")
 cat("\n")
}



plot.balance <- function(x, longcovnames=NULL,
                main="Standardized Difference in Means",
                v.axis=TRUE,
                cex.main=1, cex.vars=0.8, cex.pts=0.8,
                mar=c(0, 3, 5.1, 2), plot=TRUE, ...)
{
  K <- dim(x$diff.means.raw)[1]
  idx <- 1:K

  covnames <- x$covnames
  # prepare for plot use

  est <- x$diff.means.raw[,3]
  sd <- x$diff.means.raw[,6]
  est2 <- x$diff.means.matched[,3]
  sd2 <- x$diff.means.matched[,6]

  # x.range <- range (c(est,est2)/c(sd,sd2))
  # x.range[2] <- x.range[2] +.3
  # A <- -x.range[1]/(x.range[2]-x.range[1])
  # B <- 1/(x.range[2]-x.range[1])
  # pts <- A + B*(est/sd)              # before matched.dat
  # pts2 <- A + B*(est2/sd2)           # after macthed

  pts <-  est/sd                      # before matched.dat
  pts2 <- est2/sd2                    # after macthed
  #x.range <- c(jitter(min(c(pts, pts2)),15), max(c(pts,pts2)+.105))


  # tune the graphic console
  #par (mar=mar, mgp=mgp, oma=oma, tcl=tcl)


  par(mar = c(0, 3, 5.1, 2))
  if (is.null(longcovnames)) {
      longcovnames <- covnames
      maxchar <- max(sapply(longcovnames, nchar))
  }
  else {
      maxchar <- max(sapply(longcovnames, nchar))
  }
  min.mar <- par("mar")
  mar[2] <- min(min.mar[2], trunc(mar[2] + maxchar/10)) + mar[2] + 0.1
  par(mar = mar)

  if(plot){
     # plot the estimates
     plot(c(pts,pts2), c(idx,idx),
         bty="n", xlab="", ylab="",
         xaxt="n", yaxt="n", #xaxs="i",
         #yaxs="i",
         type="n",
         #ylim=c(max(idx)+.25, min(idx)-.25),
         #xlim=x.range,
         main=main, cex.main=cex.main,...)
     abline(v=0, lty=2)
     points(pts, idx, cex=cex.pts)          # before matched
     points(pts2, idx, pch=19, cex=cex.pts) # after matched
     if (v.axis){
         axis(3)
     }
     if (is.null(longcovnames)){
         axis(2, at=1:K, labels=covnames[1:K],
             las=2, hadj=1, lty=0, cex.axis=cex.vars)
     }
     else{
         axis(2, at=1:K, labels=longcovnames[1:K],
             las=2, hadj=1, lty=0, cex.axis=cex.vars)
     }
  }
  else{
    plot(c(pts,pts2), c(idx,idx),
      bty="n", xlab="", ylab="",
      xaxt="n", yaxt="n", #xaxs="i",
      #yaxs="i",
      type="n", axes=FALSE,
      #ylim=c(max(idx)+.25, min(idx)-.25),
      #xlim=x.range,
      main="", cex.main=cex.main,...)
  }
  return(list("raw"=pts, "matched"=pts2))
}
