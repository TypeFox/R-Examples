################################
## INTERNAL UTILITY FUNCTIONS ##
################################
##INTERNAL functions in this file:
## commonPrintST
## commonSummaryST
## internalComputeLTA
## internalPlotPredictions
## internalPlotPredictChecks
## internalFindIDplot
## internalCheckIDplot
## stCheckLoglikeIn
## stCheckType
## stCheckX
## internalFixNuggetUnobs
## checkDimInternal
## createSTmodelInternalDistance
## internalQQnormPlot
## internalScatterPlot
## internalSTmodelCreateF
## internalSummaryPredCVSTmodel

########################################################################
## Parts of S3 print that are common for STdata and STmodel, INTERNAL ##
########################################################################
commonPrintST <- function(x, name, print.type, type=NULL){
  if( print.type==1 ){
    ##general information regarding number of observations
    cat( sprintf("%s-object with:\n", name) )
    cat(sprintf("\tNo. locations: %d (observed: %d)\n",
                dim(x$covars)[1], length(unique(x$obs$ID)) ))
    cat(sprintf("\tNo. time points: %d (observed: %d)\n",
                max(dim(x$trend)[1], length(unique(x$obs$date))),
                length(unique(x$obs$date)) ))
    cat(sprintf("\tNo. obs: %d\n\n",length(x$obs$obs) ))
    ##Trends
    if( is.null(x$trend) ){
      cat("No trend specified\n")
    }else{
      if( dim(x$trend)[2]==1 ){
        cat("Only constant temporal trend,")
      }else{
        cat( sprintf("Trend with %d basis function(s):\n",
                     dim(x$trend)[2]-1))
        print( names(x$trend)[ -which(names(x$trend)=="date") ] )
      }
      cat( sprintf("with dates:\n\t%s\n",
                   paste(range(x$trend$date), collapse=" to ") ))
    }##if( is.null(x$trend) ){...}else{...}
    if( !is.null(x$old.trend) ){
      cat( "Observations have been detrended.\n")
    }
    cat("\n")
  }##if( print.type==1 )
  if( print.type==2 && !is.null(type) ){
    if( length(type)!=dim(x$covars)[1] ){
      stop( paste("length(type) =",length(type),
                  "has to equal the number fo sites =",
                  dim(x$covars)[1]) )
    }
    if( !is.factor(type) ){
      warning("Attempting to coerc 'type' to factor")
      type <- as.factor(type)
    }
    cat( "All sites:")
    print( table(type,dnn="") )
    cat( "Observed:")
    print( table(type[x$covars$ID %in% unique(x$obs$ID)]) )
    cat("\n")
    for(i in  levels(type)){
      I <- (x$obs$ID %in% x$covars$ID[type==i])
      cat( sprintf("For %s:\n",i) )
      if( sum(I)!=0 ){
        cat( sprintf("  Number of obs: %d\n", sum(I)) )
        cat( sprintf("  Dates: %s\n", paste(range(x$obs$date[I]),
                                            collapse=" to ")) )
      }else{
        cat("  No observations\n")
      }
    }##for(i in  levels(type))
  }##if( print.type==2 && !is.null(type) )
}##function commonPrintST

##########################################################################
## Parts of S3 summary that are common for STdata and STmodel, INTERNAL ##
##########################################################################
commonSummaryST <- function(object, type=NULL){
  out <- list()
  if( dim(object$obs)[1]!=0 ){
    out$obs <- summary( object$obs[, c("obs","date"), drop=FALSE])
  }
  if( !is.null(object$trend) ){
    out$trend <- summary(object$trend)
  }
  ##and possibly observations by type
  if( !is.null(type) ){
    if( length(type)!=dim(object$covars)[1] ){
      stop( paste("length(type) =",length(type),
                  "has to equal the number fo sites =",
                  dim(object$covars)[1]) )
    }
    if( !is.factor(type) ){
      warning("Attempting to coerc 'type' to factor")
      type <- as.factor(type)
    }
    out$obs.by.type <- vector("list", length(levels(type)))
    names(out$obs.by.type) <- as.character(levels(type))
    for(i in levels(type)){
      I <- (object$obs$ID %in% object$covars$ID[type==i])
      if( sum(I)!=0 ){
        out$obs.by.type[[i]] <- summary(object$obs[I, c("obs","date"), drop=FALSE])
      }
    }##for(i in  levels(type))
  }##if( !is.null(type) )

  return(out)
}##function commonSummaryST

#################################################
## Common helper functions for predict.STmodel ##
#################################################
internalComputeLTA <- function(LTA, EX, T, V=NULL, V.pred=NULL,
                               E.vec=NULL, E.vec.pred=NULL){
  if( is.null(V) ){
    ##only compute mean(EX)
    res <- matrix(NA, dim(EX)[2], length(LTA))
    rownames(res) <- colnames(EX)
  }else{
    ##also compute V( mean(EX) )
    res <- matrix(NA, dim(EX)[2]+2, length(LTA))
    rownames(res) <- c(colnames(EX),"VX","VX.pred")
    if( is.null(E.vec) ){ E.vec <- rep(1,dim(V)[1]) }
    if( is.null(E.vec.pred) ){ E.vec.pred <- E.vec }
  }
  
  for(j in 1:length(LTA)){
    Ind.LTA <- T %in% LTA[[j]]
    LTA.tmp.res <- colMeans(EX[Ind.LTA,])
    
    if( !is.null(V) ){
      LTA.tmp.res <- c(LTA.tmp.res,
                       sum(E.vec[Ind.LTA] *
                           (V[Ind.LTA,Ind.LTA,drop=FALSE] %*%
                            E.vec[Ind.LTA])) / sum(Ind.LTA)^2,
                       sum(E.vec.pred[Ind.LTA] *
                           (V.pred[Ind.LTA,Ind.LTA,drop=FALSE] %*%
                            E.vec.pred[Ind.LTA])) / sum(Ind.LTA)^2)
    }
    res[,j] <- LTA.tmp.res
  }##for(j in 1:length(LTA.tmp))
  return(res)
}##function internalComputeLTA
      

#########################################################################
## Parts of S3 plot that appends suitable ylab, xlab, main to ... args ##
#########################################################################
internalPlotFixArgs <- function(args, default=NULL, add=NULL){
  if( is.null(default) || !is.list(default) ){
    return( c(args,add) )
  }
  for(i in names(default)){
    ##does the name exist in args, if not use default value.
    if( is.null(args[[i]]) ){
      args[[i]] <- default[[i]]
    }
  }
  return( c(args,add) )
}##function internalPlotFixArgs

###########################################################################
## Parts of S3 plot that are common for predictSTmodel and predCVSTmodel ##
###########################################################################
internalPlotPredictions <- function(plot.type, ID, pred, obs, col, pch, cex,
                                    lty, lwd, p, add, transform=NULL, ...){
  ##compute the quantile
  q <- qnorm((1-p)/2, lower.tail=FALSE)
  ##ensure that lty, lwd, pch, and cex are of length==2
  if( length(lty)==1 ){ lty = c(lty,lty) }
  if( length(lwd)==1 ){ lwd = c(lwd,lwd) }
  if( length(pch)==1 ){ pch = c(pch,pch) }
  if( length(cex)==1 ){ cex = c(cex,cex) }

  if( plot.type=="obs" ){
    ##drop unobserved (or unpredicted)
    I <- !is.na(obs) & !is.na(pred$x)
    pred <- pred[I,,drop=FALSE]
    obs <- obs[I]
    ##plot as a function of sorted observations
    I <- order(obs)
    ##reorder
    obs <- obs[I]
    pred <- pred[I,,drop=FALSE]
    ##use t for x-axis below
    t <- obs
    xlab <- "observations"
  }else{
    ##dates
    t <- convertCharToDate(pred$date)
    xlab <- "time"
  }
  if( length(pred$x)==0 ){
    stop( paste("No observations for ID =", paste(ID,collapse=", ")) )
  }
  
  ##compute CI:s
  if( !is.null(transform) ){
    CI.u <- exp( pred$x.log + q*pred$sd )
    CI.l <- exp( pred$x.log - q*pred$sd )
  }else{
    CI.u <- pred$x + q*pred$sd
    CI.l <- pred$x - q*pred$sd
  }
  
  if( length(ID)==1 ){
    ID.name <- ID
  }else{
    ID.name <- paste(range(ID),collapse=" to ")
  }
  ##plot the results
  if(!add){
    args <- internalPlotFixArgs(list(...),
                                default=list(main=ID.name, xlab=xlab,
                                  ylab="predictions",
                                  ylim=range(c(CI.u, CI.l, obs, pred$x),
                                    na.rm=TRUE)),
                                add=list(x=t, y=pred$x, type="n"))
    do.call(plot, args)
  }
  ##Plot the polygon, NA if missing -> no polygon
  polygon(c(t,rev(t)), c(CI.u, rev(CI.l)),
          border=col[3], col=col[3])
  ##plot the predictions
  col.1 <- col[1]
  if( col.1=="ID" ){
    col.1 <- as.double( as.factor(pred$ID) )
  }
  if( !is.na(lty[1]) )
    lines(t, pred$x, col=col.1, lty=lty[1], lwd=lwd[1])
  if( !is.na(pch[1]) )
    points(t, pred$x, col=col.1, pch=pch[1], cex=cex[1])
  if( plot.type=="time" ){
    ##and the observations
    if( !is.na(lty[2]) )
      lines(t, obs, col=col[2], lty=lty[2], lwd=lwd[2])
    if( !is.na(pch[2]) )
      points(t, obs, col=col[2], pch=pch[2], cex=cex[2])
  }else{
    if( !is.na(lty[2]) )
      abline(0, 1, lty=lty[2], col=col[2], lwd=lwd[2])
  }
  
  ##return
  return(invisible())
}##function internalPlotPredictions


###########################################################################
## Common helper functions for S3 plot.predictSTmodel/plot.predCVSTmodel ##
###########################################################################
internalPlotPredictChecks <- function(plot.type, pred.type,
                                      pred.var, transform=NULL){
  plot.type <- match.arg(arg=plot.type, choices=c("time", "obs"))
  ##check pred.type to use
  if( is.null(transform) ){
    pred.type <- match.arg(arg=pred.type,
                           choices=c("EX", "EX.mu", "EX.mu.beta"))
  }else{
    pred.type <- match.arg(arg=pred.type,
                           choices=c("EX", "EX.pred", "EX.mu", "EX.mu.beta"))
  }
  if( !(pred.type %in% c("EX","EX.pred")) ){
    pred.var <- "DO NOT USE"
  }else{
    if( pred.var ){
      pred.var <- "VX.pred"
    }else{
      pred.var <- "VX"
    }
    if( !is.null(transform) ){
      pred.var <- paste("log.", pred.var, sep="")
      if( sum( grepl(".pred$", c(pred.type,pred.var)) )==1 ){
        warning("Using VX.pred with EX (or vice versa)")
      }
      ##also add log.EX for computations of CI
      pred.type <- c(pred.type, "log.EX")
    }
  }
  return( list(plot.type=plot.type, pred.type=pred.type, pred.var=pred.var) )
}##function internalPlotPredictChecks

internalFindIDplot <- function(ID, names){
  if( is.null(ID) ){
    ID <- 1
  }
  if( is.numeric(ID) ){
    ID <- names[ID]
  }
  if( !is.character(ID) ){
    stop("ID not character, or not inferable.")
  }
  return(ID)
}##function internalFindIDplot

internalCheckIDplot <- function(ID, y){
  if( ID=="all" || length(ID)!=1 ){
    ID.all <- TRUE
  }else{
    ID.all <- FALSE
  }
  if( ID.all && y=="time" ){
    stop("ID=all (or several locations) and y=time incompatable")
  }
  return(ID.all)
}##funciton internalCheckIDplot

#########################################################################
## Check parameter input to loglikeST, and related functions, INTERNAL ##
#########################################################################
stCheckLoglikeIn <- function(x, x.fixed, type){
  stCheckType(type)
  if( !is.null(x.fixed) ){
    if( length(x)!=sum(is.na(x.fixed)) ){
      stop("length(x) must match number of NA:s in x.fixed.")
    }
    x.fixed[ is.na(x.fixed) ] <- x
    x <- x.fixed
  }
  return(x)
}

stCheckType <- function(type){
  ##check if type is valid
  if( !(type %in% c("r","p","f")) ){
    stop("Unknown option for type, valid options are (r)eml, (p)rofile, or (f)ull.")
  }
}

############################################################
## Check parameter input to estimation and MCMC, INTERNAL ##
############################################################
stCheckX <- function(x, x.fixed, dimensions, type, object){
  ##Check if we need to truncate or expand
  if(dim(x)[1] != dimensions$nparam && dim(x)[1] != dimensions$nparam.cov){
    stop( paste("dim(x)[1] must be", dimensions$nparam,
                "or", dimensions$nparam.cov) )
  }
  ##check x.fixed
  if( !is.null(x.fixed) && (!is.vector(x.fixed) || !is.numeric(x.fixed)) ){
    stop("'x.fixed' must be a numeric vector.")
  }else if( is.null(x.fixed) ){
    ##expand
    x.fixed <- rep(NA, dim(x)[1])
  }else if(length(x.fixed) != dimensions$nparam &&
           length(x.fixed) != dimensions$nparam.cov){
    stop( paste("length(x.fixed) must be", dimensions$nparam,
                "or", dimensions$nparam.cov) )
  }
 
  if( dim(x)[1]==dimensions$nparam && type!="f"){
    ##requested REML or profile but provided full parameter
    ##starting points -> truncate
    I <- (dimensions$nparam-dimensions$nparam.cov+1):dimensions$nparam
    if( length(x.fixed)==dimensions$nparam ){
      x.fixed <- x.fixed[I]
    }
    x <- x[I,,drop=FALSE]
  }else if( dim(x)[1]==dimensions$nparam.cov && type=="f" ){
    ##requested full but provided parameters for REML or profile -> expand
    x.old <- x
    x <- matrix(NA, dimensions$nparam, dim(x)[2])
    for(i in 1:dim(x)[2]){
      ##compute alpha and gamma as cond. exp. given data
      tmp <- predict.STmodel(object, x.old[,i], only.pars=TRUE, type="p")$pars
      x[,i] <- c(tmp$gamma.E, tmp$alpha.E, x.old[,i])
    }
    if( length(x.fixed)==dimensions$nparam.cov ){
      x.fixed <- c(rep(NA,dimensions$nparam-dimensions$nparam.cov), x.fixed)
    }
  }
  ##add names to x.fixed.
  if( type!="f" ){
    names(x.fixed) <- loglikeSTnames(object, all=FALSE)
  }else{
    names(x.fixed) <- loglikeSTnames(object, all=TRUE)
  }

  ##reduce x
  x.all <- x
  x <- x[is.na(x.fixed),,drop=FALSE]
  x.all[!is.na(x.fixed),] <- matrix(x.fixed[!is.na(x.fixed)],
                                    ncol=dim(x.all)[2])
  return( list(x.all=x.all, x=x, x.fixed=x.fixed) )
}##fucntion stCheckX

###################################################
## Create nugget for all sites in STmodel adding ##
## nugget.unobs.in to unobserved locations.      ##
###################################################
internalFixNuggetUnobs <- function(nugget.unobs.in, STmodel, nugget){
  ##nugget for unobserved sites
  nugget.unobs <- matrix(NA, length(STmodel$locations$ID), 1)
  I <- match(rownames(nugget), STmodel$locations$ID)
  if( !any(is.na(I)) ){
    nugget.unobs[ I ] <- nugget
  }
  if( length(nugget.unobs.in)!=1 &&
     length(nugget.unobs.in)!=sum(is.na(nugget.unobs)) ){
    stop( paste("Needs,", 1, "or", sum(is.na(nugget.unobs)),
                "elements in nugget.unobs") )
  }
  nugget.unobs[ is.na(nugget.unobs) ] <- nugget.unobs.in
  ##add names
  rownames(nugget.unobs) <- STmodel$locations$ID
  return( nugget.unobs )
}##function internalFixNuggetUnobs

###########################################################
## Check dimensions for block matrix functions, INTERNAL ##
###########################################################
checkDimInternal <- function(X){
  dim.X <- sapply(X, dim)
  if( any(dim.X[1,1]!=dim.X[1,]) ){
    stop("all elements in X must have same number of columns.")
  }
  dim <- list(n=dim.X[1,1], m=length(X), p=dim.X[2,])
  return(dim)
}##function checkDimInternal

###########################################################
## Compute distance-matrices for createSTmodel, INTERNAL ##
###########################################################
createSTmodelInternalDistance <- function(STmodel){
  if( dim(STmodel$obs)[1]!=0 ){
    I.idx <- unique(STmodel$obs$idx)
    ##calculate distance matrices (only for observed locations),
    ##different for nu and beta
    I.idx <- 1:max(I.idx)
    STmodel$D.nu <- crossDist(STmodel$locations[I.idx, c("x.nu","y.nu"),
                                                drop=FALSE])
    colnames(STmodel$D.nu) <- rownames(STmodel$D.nu) <- STmodel$locations$ID[I.idx]
    STmodel$D.beta <- crossDist(STmodel$locations[I.idx, c("x.beta","y.beta"),
                                                  drop=FALSE])
    colnames(STmodel$D.beta) <- rownames(STmodel$D.beta) <- STmodel$locations$ID[I.idx]
    
    ##count number of observations at each time point
    dates <- sort(unique(STmodel$obs$date))
    STmodel$nt <- double(length(dates))
    for(i in c(1:length(dates))){
      STmodel$nt[i] <- sum( STmodel$obs$date==dates[i] )
    }
  }
  return( STmodel )
}##function createSTmodelInternalDistance

########################################################################
## Common helper functions for S3 qqnorm.STdata/STmodel/predCVSTmodel ##
########################################################################
internalQQnormPlot <- function(Y, ID, main, group, col, norm, line, ...){
  ##check inputs, first ID
  ID.unique <- unique(Y$ID)
  ID <- internalFindIDplot(ID, ID.unique)
  ID.all <- internalCheckIDplot(ID, "obs")

  ##if ID="all", find all possible ID names
  if(ID.all && length(ID)==1 && ID=="all"){
    ID <- ID.unique
  }
  ##extract residuals
  Ind <- Y$ID %in% ID
  Y <- Y$obs[Ind]
  
  ##extract group and colour information
  if( length(group)==length(Ind) ){
    group <- group[Ind]
  }else if( length(group)!=0 && length(group)!=sum(Ind) ){
    stop( sprintf("length(group) should be 0, %d, or %d", length(Ind), sum(Ind)) )
  }
  if( length(col)==1 ){
    col <- rep(col, sum(Ind))
  }else if( length(col)==length(Ind) ){
    col <- col[Ind]
  }else if( length(col)!=sum(Ind) ){
    stop( sprintf("length(col) should be 1, %d, or %d", length(Ind), sum(Ind)) )
  }

  ##standard QQ-plot
  qqnorm(Y, col=col, main=main, ...)
  if( norm ){
    abline(0, 1, lty=1)
  }
  if( line!=0 ){
    qqline(Y, lty=line)
  }

  ##QQ-plot for each sub-group
  if( !is.null(group) && length(unique(group))>1 ){
    for(i in unique(group)){
      qqnorm(Y[group==i], col=col[group==i], main=paste(main, i), ...)
      if( norm ){
        abline(0, 1, lty=1)
      }
      if( line!=0 ){
        qqline(Y[group==i], lty=line)
      }
    }##for(i in unique(group))
  }##if( !is.null(group) && length(unique(group))>1 )

  return(invisible())
}##function internalQQnormPlot

#############################################################################
## Common helper functions for S3 scatterPlot.STdata/STmodel/predCVSTmodel ##
#############################################################################
internalScatterPlot <- function(obs, covar, trend, data, subset, group, pch,
                                col, cex, lty, add, smooth.args, ...){
  ##need to specify either covar or trend
  if( is.null(covar) && is.null(trend) ){
    message("Both covar and trend are NULL, assuming trend=1")
    trend <- 1
  }
  if( !is.null(covar) && !is.null(trend) ){
    warning("Both covar and trend specified, ignoring trend")
    trend <- NULL
  }
  if( !is.null(group) ){
    ##ensure group is a factor
    group <- as.factor(group)
    if( length(group)!=dim(obs)[1] ){
      stop("length(group) must match the number of observations.")
    }
  }else{
    group <- as.factor( rep(1, dim(obs)[1]) )
  }#if( !is.null(group) ){...}else{...}
  ##adjust pch/col/cex/lty
  if( length(pch)==1 ){
    pch <- rep(pch, nlevels(group))
  }else if( length(pch)!=nlevels(group) ){
    stop("length(pch) should be 1 or nlevels(group)")
  }
  if( length(col)==1 ){
    col <- rep(col, nlevels(group)+1)
  }else if( length(col) != (nlevels(group)+1) ){
    stop("length(col) should be 1 or nlevels(group)+1")
  }
  if( length(cex)==1 ){
    cex <- rep(cex, nlevels(group))
  }else if( length(cex)!=nlevels(group) ){
    stop("length(cex) should be 1 or nlevels(group)")
  }
  if( length(lty)==1 ){
    lty <- rep(lty, nlevels(group)+1)
  }else if( length(lty) != (nlevels(group)+1) ){
    stop("length(lty) should be 1 or nlevels(group)+1")
  }

  ##subset the data
  if( !is.null(subset) ){
    I <- obs$ID %in% subset
    obs <- obs[I,,drop=FALSE]
    if( dim(obs)[1]==0 ){
      stop("No observations left after subsettting")
    }
    ##and the grouping
    if( !is.null(group) ){ group <- group[I] }
  }

  ##extract x and y for plotting (keeping names to get xlab/ylab right)
  y <- obs[,1,drop=FALSE]
  if( !is.null(covar) ){
    x <- data$covars[,covar,drop=FALSE]
    x <- x[match(obs$ID, data$covars$ID),,drop=FALSE]
    ##coerc x to numeric (i.e. first, and only, element of the data.frame)
    x[[1]] <- as.numeric(x[[1]])
  }else{
    if( trend==0 ){
      x <- data.frame(const=rep(1, length(y)))
    }else{
      x <- data$trend[,trend,drop=FALSE]
      x <- x[match(obs$date, data$trend$date),,drop=FALSE]
    }
  }
  ##stack x,y and create vectors fro loess.smooth
  XY <- cbind(x,y)
  x <- x[,1]
  y <- y[,1]

  ##plotting
  if(add==FALSE){
    plot(XY, type="n", ...)
  }
  points(XY, pch=pch[group], col=col[group], cex=cex[group])
  
  ##add smooths
  if( nlevels(group)!=1 ){
    for(i in 1:nlevels(group)){
      I <- group==levels(group)[i]
      if( !is.na(lty[i]) && sum(I)>0){
        prediction <- do.call(loess.smooth,
                              args=c(list(x=x[I], y=y[I]),
                                smooth.args))
        lines(prediction, col=col[i], lty=lty[i])
      }
    }
  }
  if( !is.na(lty[nlevels(group)+1]) ){
    prediction <- do.call(loess.smooth, args=c(list(x=x, y=y),
                                          smooth.args))
    lines(prediction, col=col[nlevels(group)+1], lty=lty[nlevels(group)+1])
  }
  
  return(invisible())
}##function internalScatterPlot

############################################################################
## Common helper functions createSTmodel and updateTrend.STmodel, fixes F ##
############################################################################
internalSTmodelCreateF <- function(STmodel){
  ##match times with observations
  F <- STmodel$trend[ match(STmodel$obs$date, STmodel$trend$date),,drop=FALSE]
  ##add intercept
  if( dim(F)[1]!=0 ){
    F <- cbind(1, F)
  }else{
    ##edge case for no data models (used for prediction)
    F <- cbind(matrix(0,0,1), F)
  }
  names(F)[1] <- "const"
  ##drop date column
  F <- F[,-which(names(F)=="date"),drop=FALSE]
  ##cast to matrix
  F <- data.matrix(F)
  rownames(F) <- as.character(STmodel$obs$date)
  return(F)
}##function internalSTmodelCreateF

###################################################################
## Common helper functions for predictCV, computes CV-statistics ##
###################################################################
internalSummaryPredCVSTmodel <- function(pred.struct, EX.names, transform,
                                         opts, I.n, out){
  obs <- transform( pred.struct$obs )
  EX.all <- lapply( pred.struct[EX.names], transform)
  ##and normalised residuals (these differ for transform and not transformed)
  if( !is.null(opts$transform) ){
    res <- (log(pred.struct$obs) - pred.struct$log.EX) /
      sqrt(pred.struct$log.VX.pred)
  }else if( opts$pred.var ){
    res <- pred.struct$res.norm
  }
  
  ##drop NA:s
  I <- apply(sapply(EX.all, is.na), 1, any)
  
  ##compute stats for raw observations
  for(i in EX.names){
    out$RMSE[I.n, i] <- sqrt(mean( (obs[!I] - EX.all[[i]][!I])^2 ))
  }
  out$R2[I.n, ] <- 1 - out$RMSE[I.n, ]^2 / var( obs[!I] )
  
  if( opts$pred.var || !is.null(opts$transform) ){
    ##recompute p for a two-sided CI
    q <- qnorm( (1+out$p)/2 )
    ##compute coverage
    out$coverage[I.n,1] <- mean( abs(res[!I]) < q )
  }
  return(out)
}##function internalSummaryPredCVSTmodel


