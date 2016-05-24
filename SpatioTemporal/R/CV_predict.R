###################################
## Functions for crossvalidation ##
###################################
##Functions in this file:
## predictCV.STmodel            EX:with estimateCV.STmodel
## predictCV                    EX:S3 method
## print.predCVSTmodel          EX:ok
## summary.predCVSTmodel        EX:ok
## print.summary.predCVSTmodel  EX:with summary.predCVSTmodel
## plot.predCVSTmodel           EX:with plot.predictSTmodel
## qqnorm.predCVSTmodel         EX:with qqnorm.STdata
## scatterPlot.predCVSTmodel    EX:with scatterPlot.STdata

##' @param silent Show status after each iteration?
##' @param LTA \code{TRUE}/\code{FALSE}, compute long-term temporal averages,
##'   similar to \code{\link{computeLTA}}, but with the option of including
##'   the uncertainty; see \code{\link{predict.STmodel}}.
##' @rdname estimateCV.STmodel
##' @family predCVSTmodel functions
##' @family predCVSTmodel methods
##' @method predictCV STmodel
##' @export
predictCV.STmodel <- function(object, x, Ind.cv=NULL, ..., silent=TRUE, LTA=FALSE){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")

  ##get parameters (if x is an object)
  if( inherits(x, "estCVSTmodel") ){
    if( missing(Ind.cv) || is.null(Ind.cv) ){
      Ind.cv <- x$Ind.cv
    }
    x <- coef(x, "all")
  }else if( inherits(x, "estimateSTmodel") ){
    x <- coef(x, "all")$par
  }##if( inherits(x,"estimateSTmodel") ){...}else
  ##if( inherits(x,"estimateSTmodel") ){...}
  
  ##check cross-validation groups (and force to matrix if possible)
  Ind.cv <- stCheckInternalCV(Ind.cv)
  
  ##ensure that Ind.cv is a matrix
  Ind.cv <- as.matrix(Ind.cv)
  if( dim(Ind.cv)[2]==1 ){
    N.CV.sets <- max(Ind.cv,na.rm=TRUE)
  }else{
    stop("Some observation(s) are left out in several Cv-groups.")
  }
                         
  ##expand the parameter vector so that we have one
  ##vector for each of the cv-sets we want to try.
  x <- as.matrix(x)
  if(dim(x)[2]==1){
    x <- matrix(x, length(x), N.CV.sets)
  }else if(dim(x)[2] != N.CV.sets){
    stop("Number of parameters does not match the number of cv-sets.")
  }

  ##loop over all locations and predict
  pred <- list()
  for(i in 1:N.CV.sets){
    if(!silent)
      message( sprintf("Predicting cv-set %d/%d", i, N.CV.sets) )
    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==i
    }else{
      Ind.current <- as.logical( Ind.cv[,i] )
    }
    ##create data matrices that contains observed
    object.obs <- dropObservations(object, Ind.current)
    ##and locations to be predicted
    suppressWarnings( object.pred <- dropObservations(object, !Ind.current) )

    ##compute nugget for unobserved
    nugget.unobs <- loglikeSTgetPars(x[,i], object)$cov.nu$nugget
    nugget.unobs <- nugget.unobs[object.pred$locations$ID, , drop=FALSE]
    ##should we compute LTA for the observations?
    if( LTA ){
      LTA.pred <- with(object.pred$obs, split(date, ID))
    }else{
      LTA.pred <- FALSE
    }
    ##lets obtain the predection for this set
    pred[[i]] <- predict(object.obs, x[,i], STdata=object.pred,
                         nugget.unobs=nugget.unobs, only.pars=FALSE,
                         combine.data=FALSE, LTA=LTA.pred, ...)
  }##for(i in 1:N.CV.sets){
  
  ##list containing results
  out <- list()
  #pred.obs=pred.obs, pred.by.cv=pred.by.cv, pred.all=pred.all.res)
  class(out) <- "predCVSTmodel"
  out$opts <- pred[[1]]$opts
  if( !is.null(out$opts$nugget.unobs) ){
    out$opts$nugget.unobs <- unlist(sapply(pred, function(x){
      x$opts$nugget.unobs}))
    names(out$opts$nugget.unobs) <- unlist(sapply(pred,function(x){
      rownames(x$opts$nugget.unobs)}))
  }
  ##fix the  LTA.list
  if( out$opts$LTA ){
    out$opts$LTA.list <- unlist(lapply(pred, function(x){ x$opts$LTA.list }),
                                recursive=FALSE)
  }else{
    out$opts$LTA.list <- NULL
  }
##save the CV-division
  out$Ind.cv <- Ind.cv

  ##collect results.
  ##construct a matrix matching mesa.data$obs$obs with the obs
  out$pred.obs <- object$obs[,c("obs","date","ID"), drop=FALSE]
  if( !is.null(out$opts$transform) ){
    out$pred.obs$obs <- exp( out$pred.obs$obs )
  }
  ##add some more columns
  EX.names <- names(pred[[1]])[grep("EX",names(pred[[1]]))]
  for(i in EX.names){
    out$pred.obs[[i]] <- NA
  }
  ##do we have prediction variances?
  if( out$opts$pred.var || !is.null(out$opts$transform) ){
    VX.names <- names(pred[[1]])[grep("VX|MSPE",names(pred[[1]]))]
    for(i in VX.names){
      out$pred.obs[[i]] <- NA
    }
  }else{
    VX.names <- NULL
  }
  ##loop over the CV-sets
  for(i in 1:N.CV.sets){
    Ind.current <- Ind.cv==i
    I <- pred[[i]]$I$I
    out$pred.obs[Ind.current, EX.names] <-
      sapply(pred[[i]][EX.names], function(x){x[I]})
    if( !is.null(VX.names) ){
      out$pred.obs[Ind.current, VX.names] <-
        sapply(pred[[i]][VX.names], function(x){x[I]})
    }
  }
  ##residuals
  out$pred.obs$res <- out$pred.obs$obs - out$pred.obs$EX
  ##normalised residuals
  if( out$opts$pred.var && is.null(out$opts$transform) ){
    out$pred.obs$res.norm <- out$pred.obs$res / sqrt(out$pred.obs$VX.pred)
  }
  ##LTA
  if( out$opts$LTA ){
    out$pred.LTA <- do.call(rbind,lapply(pred, function(x){x$LTA}))
    LTA.mean <- sapply(with(out$pred.obs, split(obs,ID)), mean)
    out$pred.LTA <- cbind(LTA.mean[match(rownames(out$pred.LTA),
                                         names(LTA.mean))],
                          rownames(out$pred.LTA),
                          out$pred.LTA)
    names(out$pred.LTA)[1:2] <- c("obs","ID")
    out$pred.LTA$ID <- as.character(out$pred.LTA$ID)
    rownames(out$pred.LTA) <- NULL

    ##residuals
    if( "EX.pred" %in% names(out$pred.LTA) ){
      out$pred.LTA$res <- out$pred.LTA$obs - out$pred.LTA$EX.pred
    }else{
      out$pred.LTA$res <- out$pred.LTA$obs - out$pred.LTA$EX
    }
    ##normalised residuals
    if( out$opts$pred.var ){
      out$pred.LTA$res.norm <- out$pred.LTA$res / sqrt(out$pred.LTA$VX.pred)
    }
  }##if( out$opts$LTA )
  
  ##determine if the same location is in several cross-validation groups
  if( out$opts$only.obs ){
    any.duplicates <- any( duplicated( unlist( sapply(pred, function(x){
      unique(x$I$ID)}) ) ) )
  }else{
    any.duplicates <- any( duplicated( unlist( sapply(pred, function(x){
      colnames(x$EX)}) ) ) )
  }

  ##Extract by location predictions
  if( !any.duplicates ){
    out$pred.all <- list()

    ##matrices that contain predictions at all times and locations
    out$pred.all$EX.mu <- matrix(NA, dim(object$trend)[1],
                                 dim(object$locations)[1])
    colnames(out$pred.all$EX.mu) <- object$locations$ID
    rownames(out$pred.all$EX.mu) <- as.character(object$trend$date)
    for(i in EX.names[EX.names!="EX.mu"]){
      out$pred.all[[i]] <- out$pred.all$EX.mu
    }
    if( !is.null(VX.names) ){
      for(i in VX.names){
        out$pred.all[[i]] <- out$pred.all$EX.mu
      }
    }

    ##matrices that contain beta predictions
    out$pred.all$beta <- list()
    out$pred.all$beta$mu <- matrix(NA, dim(object$locations)[1],
                                length(object$LUR))
    colnames(out$pred.all$beta$mu) <- names(object$LUR)
    rownames(out$pred.all$beta$mu) <- object$locations$ID
    out$pred.all$beta$EX <- out$pred.all$beta$mu
    if( out$opts$pred.var ){
      out$pred.all$beta$VX <- out$pred.all$beta$EX
    }

    for(i in 1:N.CV.sets){
      ##sites predicted in this CV-group
      if( out$opts$only.obs ){
        ID.names <- rownames(pred[[i]]$beta$EX)
        I <- ((match( pred[[i]]$I$ID, object$locations$ID)-1) *
              dim(out$pred.all$EX)[1] +
              match( as.character(pred[[i]]$I$date),
                    rownames(out$pred.all$EX)) )
      }else{
        ID.names <- colnames(pred[[i]]$EX)
        I <- (rep(match(ID.names, object$locations$ID)-1,
                  each=dim(out$pred.all$EX)[1]) * dim(out$pred.all$EX)[1] +
              rep(match(rownames(pred[[i]]$EX), rownames(out$pred.all$EX)),
                  length(ID.names)) )
      }
      ##extract predictions
      for(j in EX.names){
        out$pred.all[[j]][I] <- pred[[i]][[j]]
      }
      ##...variances
      if( !is.null(VX.names) ){
        for(j in VX.names){
          out$pred.all[[j]][I] <- pred[[i]][[j]]
        }
      }
      ##...and beta-fields
      out$pred.all$beta$mu[ID.names,] <- pred[[i]]$beta$mu
      out$pred.all$beta$EX[ID.names,] <- pred[[i]]$beta$EX
      if( out$opts$pred.var ){
        out$pred.all$beta$VX[ID.names,] <- pred[[i]]$beta$VX
      }
    }##for(i in 1:N.CV.sets)
  }else{
    out$pred.all.by.cv <- vector("list", N.CV.sets)
    for(i in 1:N.CV.sets){
      if( out$opts$only.obs ){
        out$pred.all.by.cv[[i]] <- list()
        ##create data matrix for the different parts
        EX.tmp <- vector("list", length(EX.names))
        names(EX.tmp) <- EX.names
        for(j in EX.names){
          EX.tmp[[j]] <- createDataMatrix(obs=pred[[i]][[j]],
                                          date=pred[[i]]$I$date, 
                                          ID=pred[[i]]$I$ID)
        }
        if( !is.null(VX.names) ){
          VX.tmp <- vector("list", length(VX.names))
          names(VX.tmp) <- VX.names
          for(j in VX.names){
            VX.tmp[[j]] <- createDataMatrix(obs=pred[[i]][[j]],
                                            date=pred[[i]]$I$date,
                                            ID=pred[[i]]$I$ID)
          }
        }
        
        ##create target matrices with ALL dates
        out$pred.all.by.cv[[i]][[ EX.names[1] ]] <-
          matrix(NA, dim(object$trend)[1], dim(EX.tmp[[1]])[2])
        colnames(out$pred.all.by.cv[[i]][[ EX.names[1] ]]) <-
          colnames(EX.tmp[[1]])
        rownames(out$pred.all.by.cv[[i]][[ EX.names[1] ]]) <-
          as.character(object$trend$date)
        for(j in EX.names[-1]){
          out$pred.all.by.cv[[i]][[j]] <- out$pred.all.by.cv[[i]][[ EX.names[1] ]]
        }
        if( !is.null(VX.names) ){
          for(j in VX.names){
            out$pred.all.by.cv[[i]][[j]] <- out$pred.all.by.cv[[i]][[ EX.names[1] ]]
          }
        }
        
        ##fit predictions into matrices
        for(j in EX.names){
          out$pred.all.by.cv[[i]][[j]][rownames(EX.tmp[[j]]),] <- EX.tmp[[j]]
        }
        if( !is.null(VX.names) ){
          for(j in VX.names){
            out$pred.all.by.cv[[i]][[j]][rownames(VX.tmp[[j]]),] <- VX.tmp[[j]]
          }
        }
        
        ##just copy beta predictions
        out$pred.all.by.cv[[i]]$beta <- pred[[i]]$beta
      }else{
        ##just copy
        pick.names <- c(EX.names, VX.names, "beta")
        out$pred.all.by.cv[[i]] <- pred[[i]][pick.names]
      }##if( out$opts$only.obs ){...}else{...}
    }##for(i in 1:N.CV.sets)
  }##if( !any(duplicated(...)) ){...}else{...}

  ##return results
  return( out )
}##function predictCV.STmodel

######################################
## General S3 methods for predictCV ##
######################################

##' @rdname estimateCV.STmodel
##' @export
predictCV <- function(object, x, Ind.cv, ...){
  UseMethod("predictCV")
}

##################################
## S3-METHODS FOR predCVSTmodel ##
##################################

##' \code{\link[base:print]{print}} method for class \code{predCVSTmodel}.
##'
##' @title Print details for \code{predCVSTmodel} object
##' @param x \code{predCVSTmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @examples
##' ##load some data
##' data(pred.cv.mesa)
##' ##print basic information for the CV-predictions
##' print(pred.cv.mesa)
##' 
##' @author Johan Lindström
##' 
##' @family predCVSTmodel methods
##' @method print predCVSTmodel
##' @export
print.predCVSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "predCVSTmodel", name="x")

  N.cv <- max(x$Ind.cv)
  by.CV <- !is.null(x$pred.all.by.cv)
  N.loc <- dim(x$pred.all$EX)
  
  cat("Cross-validation prediction for STmodel with", N.cv, "CV-groups.\n")
  if( x$opts$only.obs ){
    cat("  Predictions for observations only.\n")
  }else{
    cat("  Predictions for", N.loc[2], "locations and",
        N.loc[1], "time points.\n")
  }
  if( !is.null(x$opts$transform) ){
    cat("\tlog-Gaussian field of type:",
        x$opts$transform, "\n")
  }
  if( isTRUE(x$opts$LTA)  ){
    cat("  Temporal averages predicted for", dim(x$pred.LTA)[1],
        "locations.\n")
  }
  if( x$opts$pred.var ){
    cat("  Variances have been computed.\n")
    if( x$opts$type=="r" ){
      cat("\tRegression parameters are assumed to be unknown,\n")
      cat("\tprediction variances include uncertainties\n")
      cat("\tin regression parameters.\n\n")
    }else{
      cat("\tRegression parameters are assumed to be known,\n")
      cat("\tprediction variances do NOT include\n")
      cat("\tuncertainties in regression parameters.\n\n")
    }
  }else{
    cat("  Variances have NOT been computed.\n")
  }
  if( by.CV ){
    cat("  Same location occurs in several CV-groups.\n")
  }

  return(invisible())
}##function print.predCVSTmodel

##' \code{\link[base:summary]{summary}} method for class \code{predCVSTmodel}.
##'
##' Computes summary statistics for cross validation. Statistics that are
##' computed include RMSE, R2, and coverage of CI:s; both for all observations
##' and (possibly) stratified by date.
##' 
##' @title Computes summary details for \code{predCVSTmodel} object
##' 
##' @param object \code{predCVSTmodel} object to compute summary information
##'   for; the output from \code{\link{predictCV.STmodel}}.
##' @param pred.naive Result of naive prediction; used to compute
##'   modified R2 values. The output from \code{\link{predictNaive}}.
##' @param by.date Compute individual cross-validation statistics for each
##'   time-point. May lead to \emph{very many} statistics.
##' @param p Approximate coverage of the computed confidence bands; the
##'   confidence bands are used when computing coverage of the
##'   cross-validated predictions.
##' @param transform Transform observations and predictions (\emph{without} bias
##'   correction) \emph{before}
##'   computing statistics; see also \code{\link{computeLTA}}. \emph{Redundant}
##'   if option \code{transform} was used in \code{\link{predictCV.STmodel}} (as
##'   pass through argument to \code{\link{predict.STmodel}})
##' @param LTA Compute cross-validation statistics for the long term averages at
##'   each site, uses \code{\link{computeLTA}} to compute the averages.
##'   \code{transform} is passed to \code{\link{computeLTA}}. This is
##'   \emph{redundant} if option \code{LTA=TRUE} was uses in
##'   \code{\link{predictCV.STmodel}}.
##' @param ... Ignored additional arguments.
##' 
##' @return A \code{summary.predCVSTmodel} object.
##'
##' @examples
##' ##load some data
##' data(pred.cv.mesa)
##' 
##' ##basic summary statistics
##' summary(pred.cv.mesa)
##' 
##' @author Johan Lindström
##' 
##' @family predCVSTmodel methods
##' @method summary predCVSTmodel
##' @export
summary.predCVSTmodel <- function(object, pred.naive=NULL, by.date=FALSE,
                                  p=0.95, transform=function(x){return(x)},
                                  LTA=FALSE, ...){
  ##check class belonging
  stCheckClass(object, "predCVSTmodel", name="object")
  ##inputs
  if(p<=0 || p>=1){
    stop("'p' not a probability.")
  }
  if( !is.function(transform) ){
    stop("'transform' should be a function")
  }
  ##backwards compability with old predCVSTmodel objects
  object$opts$LTA <- isTRUE(object$opts$LTA)
  if( object$opts$LTA || !is.null(object$opts$transform) ){
    if( LTA || !missing(transform) ){
      warning("LTA and/or transform already computed by predictCV.\nOptions LTA and transform have been ignored")
    }
    LTA <- FALSE
    transform <- function(x){return(x)}
  }
  compute.LTA <- object$opts$LTA || LTA

  ##all EX-computed
  EX.names <- names(object$pred.obs)[ grep("^EX",names(object$pred.obs)) ]
  
  ##compute predictions by date?
  if( by.date ){
    pred.by.date <- split(object$pred.obs, object$pred.obs$date)
    ##and keep only dates with some predictions
    pred.by.date <- pred.by.date[!sapply(pred.by.date, function(x){
      all(is.na(x$EX)) })]
  }else{
    pred.by.date <- NULL
  }##if( by.date ){...}else{...}
  
  ##total number of summaryStats to compute
  Nstat.naive <- Nstat <- 1 + compute.LTA + length(pred.by.date)
  if( !is.null(pred.naive) ){
     Nstat.naive <- Nstat.naive + (dim(pred.naive$pred)[2]-2)
   }

  ##return the computed stats
  out <- list(RMSE=NA, R2=NA, p=p, stats=NA, coverage=NA)
  class(out) <- "summary.predCVSTmodel"
  
  ##allocate memory for each of these
  out$R2 <- matrix(NA, Nstat.naive, length(EX.names))
  out$RMSE <- matrix(NA, Nstat, length(EX.names))
  out$coverage <- matrix(NA, Nstat, 1)
  ##add names to matrix
  colnames(out$coverage) <- switch(is.null(object$opts$transform)+1,"EX.pred","EX")
  colnames(out$RMSE) <- colnames(out$R2) <- EX.names
  Rn <- "obs"
  if(compute.LTA){
    Rn <- c(Rn,"average")
  }
  Rn <- c(Rn, names(pred.by.date))
  rownames(out$RMSE) <- rownames(out$coverage) <- Rn
  ##also add naive predictions for R2
  if( !is.null(pred.naive) ){
    Rn <- c(Rn, names(pred.naive$pred)[!(names(pred.naive$pred) %in%
                                         c("ID","date"))])
  }
  rownames(out$R2) <- Rn

  ##compute RMSE and R2
  out <- internalSummaryPredCVSTmodel(object$pred.obs, EX.names, transform,
                                      opts=object$opts, I.n="obs", out)

  ##compute stats for long term average
  if( compute.LTA ){
    ##options for LTA case
    opts.tmp <- object$opts
    if( object$opts$LTA ){
      lta.tmp <- object$pred.LTA
    }else{
      lta.tmp <- computeLTA(object, transform)
      ##variance NOT computed by computeLTA
      opts.tmp$pred.var <- FALSE
    }
    ##drop transform, handed above
    opts.tmp$transform <- NULL
    
    out <- internalSummaryPredCVSTmodel(lta.tmp, EX.names,
                                        transform=function(x){return(x)},
                                        opts=opts.tmp, I.n="average", out)
    N.LTA <- dim( lta.tmp )[1]
  }##if( compute.LTA )
  
  ##compute stats for each date
  if( length(pred.by.date)>0 ){
    for( i in 1:length(pred.by.date) ){
      I.n <- names(pred.by.date)[i]
      out <- internalSummaryPredCVSTmodel(pred.by.date[[i]], EX.names,
                                          transform, object$opts, I.n, out)
    }## for( i in 1:length(pred.by.date) )
  }## if( length(pred.by.date)>0 )

  ##compute modified R2 for the naive predictions
  if( !is.null(pred.naive) ){
    ##check that pred.naive matches my predictions
    if( any( object$pred.obs$date != pred.naive$pred$date ) ){
      stop("Missmatch between dates in 'object' and 'pred.naive'")
    }
    if( any( object$pred.obs$ID != pred.naive$pred$ID ) ){
      stop("Missmatch between ID:s in 'object' and 'pred.naive'")
    }

    ##extract naive predictions
    I.naive <- !(names(pred.naive$pred) %in% c("ID","date"))
    if( is.null(object$opts$transform) ){
      pred <- transform( pred.naive$pred[,I.naive,drop=FALSE] )
      ##extract observations
      obs <- transform( object$pred.obs$obs )
    }else{
      ##fixed exp-transform
      pred <- exp( pred.naive$pred[,I.naive,drop=FALSE] )
      ##extract observations
      obs <- object$pred.obs$obs
    }

    ##check for which points we have predictions
    ##(same condition as above EX should exist)
    I <- !apply(is.na(object$pred.obs[,EX.names]), 1, any)
    ##compute modified R2:s
    for(i in 1:dim(pred)[2]){
      V.naive <- mean( (pred[I,i] - obs[I])^2 )
      out$R2[names(pred)[i], ] <- 1 - out$RMSE["obs",]^2 / V.naive
    }
  }
  
  ##ensure that R2 is >= 0
  out$R2 <- pmax(out$R2, 0)

  ##some CV-statistics
  out$stats <- list(N.cv=max(object$Ind.cv), Npred=sum(!is.na(object$pred.obs$EX)),
                    pred.var=object$opts$pred.var)
  if( compute.LTA ){
    out$stats$N.LTA <- N.LTA
  }
  if( !object$opts$pred.var && is.null(object$opts$transform) ){
    out$coverage <- NULL
  }
  out$stats$pred.var <- !is.null(out$coverage)
  
  return(out)
}##function summary.predCVSTmodel

##' \code{\link[base:print]{print}} method for class \code{summary.predCVSTmodel}.
##'
##' @title Print details for \code{summary.predCVSTmodel} object
##' @param x \code{summary.predCVSTmodel} object to print information for.
##' @param ... Additional arguments, passed to
##'   \code{\link[base:print]{print.table}}.
##' @return Nothing
##'
##' @author Johan Lindström
##' 
##' @family predCVSTmodel methods
##' @method print summary.predCVSTmodel
##' @export
print.summary.predCVSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "summary.predCVSTmodel", name="x")

  cat("Cross-validation predictions for STmodel with", x$stats$N.cv,
      "CV-groups.\n")
  cat("  Predictions for", x$stats$Npred, "observations.\n")
  if( !is.null(x$stats$N.LTA) ){
    cat("  Temporal averages for", x$stats$N.LTA, "locations.\n\n")
  }

  cat("RMSE:\n")
  print( x$RMSE )
  cat("\n")
  
  cat("R2:\n")
  print( x$R2 )
  cat("\n")

  if( x$stats$pred.var ){
    cat("Coverage of ", 100*x$p, "% prediction intervals:\n", sep="")
    print( x$coverage )
    cat("\n")
  }else{
    cat("Coverage (i.e. prediciton variance) have NOT been computed.\n")
  }
  cat("\n")

  return(invisible())
}##function print.summary.predCVSTmodel

##' @rdname plot.predictSTmodel
##' @family predCVSTmodel methods
##' @importFrom graphics plot
##' @method plot predCVSTmodel
##' @export
plot.predCVSTmodel <- function(x, y="time", ID=colnames(x$pred.all$EX)[1],
                                col=c("black","red","grey"), pch=c(NA,NA),
                                cex=c(1,1), lty=c(1,1), lwd=c(1,1), p=0.95,
                                pred.type="EX", pred.var=TRUE,
                                add=FALSE, ...){
  ##check class belonging
  stCheckClass(x, "predCVSTmodel", name="x")
  ##we have to use y, cast to resonable name
  tmp <- internalPlotPredictChecks(y, pred.type, pred.var, x$opts$transform)
  plot.type <- tmp$plot.type
  pred.type <- tmp$pred.type
  pred.var <- tmp$pred.var
  
  ##check ID
  ID <- internalFindIDplot(ID, colnames(x$pred.all$EX))
  ID.all <- internalCheckIDplot(ID, y)

  ##plotting options
  if( !ID.all ){
    ##pick out predictions, two cases
    date <- rownames(x$pred.all[[pred.type[1]]])
    sd <- x$pred.all[[pred.var]][,ID]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x$pred.all[[ pred.type[1]] ][,ID], sd=sd,
                       date=date, x.log=NA, stringsAsFactors=FALSE)
    if( length(pred.type)==2 ){
      pred$x.log <- x$pred.all[[ pred.type[2] ]][,ID]
    }
    ##pick out observations, if any
    obs <- x$pred.obs[x$pred.obs$ID==ID,,drop=FALSE]
    if( !is.null(obs) ){
      tmp <- obs
      obs <- matrix(NA, dim(pred)[1], 1)
      rownames(obs) <- as.character(pred$date)
      obs[as.character(tmp$date),] <- tmp$obs
    }
  }else{
    ##all only relevant when we have plot by observations
    ##pick out predictions, two cases
    if(length(ID)==1 && ID=="all"){
      ID <- unique(x$pred.obs$ID)
    }
    I.ID <- x$pred.obs$ID %in% ID
    sd <- x$pred.obs[I.ID,pred.var]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x$pred.obs[I.ID,pred.type[1]], sd=sd,
                       date=x$pred.obs[I.ID,"date"],
                       ID=x$pred.obs[I.ID,"ID"],
                       x.log=NA, stringsAsFactors=FALSE)
    if( length(pred.type)==2 ){
      pred$x.log <- x$pred.obs[I.ID,pred.type[2]]
    }
    ##pick out observations, if any
    obs <- x$pred.obs[I.ID,"obs"]
    
    if( dim(pred)[1]!=length(obs) ){
      stop("Number of observed locations do not match no. predicted.")
    }
  }##if( !ID.all ){...}else{...}

  internalPlotPredictions(plot.type, ID, pred, obs, col, pch, cex,
                          lty, lwd, p, add, x$opts$transform, ...)
  
  return(invisible())
}##function plot.predCVSTmodel


##' @param norm \code{TRUE}/\code{FALSE}, plot normalised (mean=0, sd=1) or raw
##'   cross-validation residuals. If \code{norm=TRUE} a 0-1 line is added, to
##'   indicate what normalised residuals should look like.
##' @param org.scale \code{TRUE}/\code{FALSE} scatter plots on the original
##'   untransformed scale, or using \code{exp(y)}. Only relevant if \code{x} was
##'   computed using \code{transform} in \code{\link{predictCV.STmodel}} (as
##'   pass through argument to \code{\link{predict.STmodel}})
##' @rdname qqnorm.STdata
##' @family predCVSTmodel methods
##' @importFrom stats qqnorm
##' @method qqnorm predCVSTmodel
##' @export
qqnorm.predCVSTmodel <- function(y, ID="all", main="Q-Q plot for CV residuals",
                                 group=NULL, col=1, norm=FALSE, line=0,
                                 org.scale=TRUE, ...){ 
  ##check class belonging
  stCheckClass(y, "predCVSTmodel", name="y")

  if( !is.null(y$opts$transform) ){
    if( org.scale ){
      ##inverse tansform for reasonable results
      y$pred.obs$res <- log(y$pred.obs$obs) - y$pred.obs$log.EX
      y$pred.obs$res.norm <- y$pred.obs$res / sqrt(y$pred.obs$log.VX.pred)
    }else{
      warning("log-Gaussian field: Residuals are unlikely to be Gaussian.")
    }
  }
  if(norm){
    if( inherits(try(Y <- y$pred.obs[, c("ID","res.norm")], silent=TRUE),
                 "try-error") ){
      stop("Unable to use norm=TRUE; most likely variance NOT computed in predictCV")
    }
  }else{
    Y <- y$pred.obs[, c("ID","res")]
  }
  colnames(Y) <- c("ID","obs")
  
  internalQQnormPlot(Y, ID, main, group, col, norm, line, ...)
  
  return(invisible())
}##function qqnorm.predCVSTmodel

##' @param type What to use in the scatter plot, valid options are \code{"obs"}
##'   for observations, \code{"res"} residuals, and \code{"res.norm"} for
##'   normalised residuals.
##' @param STdata \code{STdata} or \code{STmodel} containing covariates and
##'   trend against which to plot.
##' @param org.scale \code{TRUE}/\code{FALSE} scatter plots on the original
##'   untransformed scale, or using \code{exp(y)}. Only relevant if \code{x} was
##'   computed using \code{transform} in \code{\link{predictCV.STmodel}} (as
##'   pass through argument to \code{\link{predict.STmodel}})
##' @rdname scatterPlot.STdata
##' @family predCVSTmodel methods
##' @method scatterPlot predCVSTmodel
##' @export
scatterPlot.predCVSTmodel <- function(x, covar=NULL, trend=NULL, pch=1, col=1,
                                      cex=1, lty=1, subset=NULL, group=NULL,
                                      add=FALSE, smooth.args=NULL, STdata,
                                      type=c("obs","res","res.norm"),
                                      org.scale=TRUE, ...){
  ##check class belonging
  stCheckClass(x, "predCVSTmodel", name="x")

  ##check type
  type <- match.arg(type)

  if( !is.null(x$opts$transform) ){
    if( org.scale ){
      ##inverse tansform for reasonable results
      x$pred.obs$obs <- log(x$pred.obs$obs)
      x$pred.obs$res <- x$pred.obs$obs - x$pred.obs$log.EX
      x$pred.obs$res.norm <- x$pred.obs$res / sqrt(x$pred.obs$log.VX.pred)
    }else{
      warning("Using log-Gaussian field, scatter plot may be strange")
    }
  }

  ##check aux.data
  stCheckClass(STdata, c("STdata","STmodel"), name="STdata")
  
  ##extract observations/residuals/norm residuals
  if( inherits(try(obs <- x$pred.obs[, c(type,"ID","date")], silent=TRUE),
               "try-error") ){
    stop("Unable to use type=", type,
         "; most likely variance NOT computed in predictCV")
  }
  
  ##extract aux.data
  if( inherits(STdata, "STdata") ){
    covars <- STdata$covars
    ##add trend to avoid future problems
    if( is.null(STdata$trend) ){
      STdata <- updateTrend(STdata, n.basis=0)
    }
  }else{
    ##if STmodel, collect location & LUR information
    covars <- cbind(STdata$locations, do.call(cbind, STdata$LUR))
    covars <- covars[,unique(colnames(covars))]
  }

  ##pass data to internalScatterPlot function
  internalScatterPlot(obs=obs, covar=covar, trend=trend, subset=subset,
                      data=list(covars=covars, trend=STdata$trend),
                      group=group, pch=pch, col=col, cex=cex,
                      lty=lty, add=add, smooth.args=smooth.args, ...)
}##function scatterPlot.predCVSTmodel
