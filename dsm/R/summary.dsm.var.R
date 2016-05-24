#' Summarize the variance of a density surface model
#'
#' Gives a brief summary of a fitted \code{dsm} variance object. 
#'
#' @aliases summary.dsm.var
#'
#' @param object a \code{dsm.var} object
#' @param alpha alpha level for confidence intervals
#' @param boxplot.coef the value of \code{coef} used to calculate the outliers
#'        see \code{\link{boxplot}}.
#' @param bootstrap.subregions list of vectors of logicals or indices for 
#'        subregions for which variances need to be calculated (only for
#'        bootstraps (see \code{\link{dsm.var.prop}} for how to use subregions
#'        with variance propagation).
#' @param \dots unused arguments for S3 compatibility
#' @return a summary object
#' @export
#'
#' @seealso dsm.var.movblk dsm.var.prop
#' @author David L. Miller
#'
summary.dsm.var<-function(object, alpha=0.05, boxplot.coef=1.5,
                  bootstrap.subregions=NULL,...){

  # storage
  sinfo<-list()

  if(object$bootstrap){
    # grab the predicted values
    #mod1.pred <- dsm.predict(object$dsm.object,
    #                         newdata=object$pred.data,
    #                         off.set=object$off.set)
    #sinfo$pred.est <- sum(mod1.pred,na.rm=TRUE)
    sinfo$pred.est <- object$study.area.total[1]

    sinfo$block.size <- object$block.size
    sinfo$n.boot <- object$n.boot
    sinfo$bootstrap <- TRUE
    sinfo$ds.uncertainty <- object$ds.uncertainty

    # bootstrap abundances
    bootstrap.abund <- object$study.area.total

    # when we don't need to do the delta method
    if(is.null(object$dsm.object$ddf) | object$ds.uncertainty){

      # if we used detection function uncertainty
      # or there was no detection function

      # variance of the bootstrap abundances is the variance
      trimmed.variance <- trim.var(bootstrap.abund[is.finite(bootstrap.abund)],
                                   boxplot.coef=boxplot.coef)

      sinfo$var <- trimmed.variance
      sinfo$se <- sqrt(trimmed.variance)

      sinfo$cv <- sinfo$se/sinfo$pred.est
      # in this case bootstrap and regular CV are the same
      sinfo$bootstrap.cv <- sinfo$cv
    }else{
      # delta method, if necessary
      #  - if we didn't incorporate detection function uncertainty in bootstrap
      #  - if there was a detection function

      ddf.summary <- summary(object$dsm.object$ddf)

      # average p standard error
      sinfo$average.p.se <- ddf.summary$average.p.se

      ## calculate the variance via the delta method
      # find the cv squared of the p
      cvp.sq <- (ddf.summary$average.p.se/
                 ddf.summary$average.p)^2
      # save that
      sinfo$detfct.cv <- sqrt(cvp.sq)

      # save the s.e. of N from bootstrap
      trimmed.variance <- trim.var(bootstrap.abund[is.finite(bootstrap.abund)],
                                   boxplot.coef=boxplot.coef)
      sinfo$N.bs.se <- sqrt(trimmed.variance)

      # cv squared of the Ns from the bootstrap
      cvNbs.sq <- (sinfo$N.bs.se/sinfo$pred.est)^2
      # save that
      sinfo$bootstrap.cv <- sqrt(cvNbs.sq)

      # cv of N
      cvN <- sqrt(cvp.sq+cvNbs.sq)
      sinfo$cv <- cvN

      # variance (delta method)
      sinfo$var <- (cvN*sinfo$pred.est)^2
      sinfo$se <- sqrt(sinfo$var)
    }

    ### general bootstrap stuff

    # how many duds did we have?
    sinfo$boxplot.coef <- boxplot.coef
    sinfo$trim.prop <- attr(trimmed.variance, "trim.prop")
    sinfo$trim.ind <- attr(trimmed.variance, "trim.ind")
    sinfo$boot.outliers <- attr(trimmed.variance, "outliers")
    sinfo$boot.infinite <- sum(is.infinite(bootstrap.abund))
    sinfo$boot.finite <- sum(!is.infinite(bootstrap.abund))
    sinfo$boot.NA <- sum(is.na(bootstrap.abund))
    sinfo$boot.NaN <- sum(is.nan(bootstrap.abund))
    sinfo$boot.usable <- sinfo$boot.finite - (sinfo$boot.outliers + 
                         sinfo$boot.infinite + sinfo$boot.NA + sinfo$boot.NaN)

    # grab the %ile c.i.s at alpha, 1-alpha and also median
    sinfo$quantiles <- quantile(bootstrap.abund[sinfo$trim.ind], 
                                c(alpha, 0.5, 1-alpha),na.rm=TRUE)
    attr(sinfo$quantiles,"names")[2] <- "Median"


    ### subregions...
    if(!is.null(bootstrap.subregions)){
      # to do the subregions, we just recurse back into this routine
      # each time we just restrict the data to those specified by
      # the corresponding entry in bootstrap.subregions

      subregions <- list()
      i<-1

      for(region.ind in bootstrap.subregions){
        # setup the subregion object
        this.object <- object
        this.object$short.var <- NULL
        this.object$study.area.total <- object$study.area.total[region.ind]
        this.object$pred.data <- object$pred.data[region.ind,]

        # summarise and store region i
        subregions[[i]] <- summary(this.object)
        i<-i+1
      }
      sinfo$subregions<-subregions
    }

  }else{
  ### varprop and "Bayesian" stuff
    sinfo$varprop <- !is.null(object$deriv)
    sinfo$saved<-object
    sinfo$bootstrap <- object$bootstrap

    # what if we had multiple areas (ie this is from a CV plot?)
    if(all(dim(as.matrix(object$pred.var))==1)){
      sinfo$se <- sqrt(object$pred.var)
    }else{
      # re run the variance calculation, putting everything together
      pd<-c()
      off<-c()
      for(i in 1:length(object$pred.data)){
        pd<-rbind(pd,object$pred.data[[i]])
        off<-rbind(off,object$off.set[[i]])
      }
      object$pred.data <- pd
      object$off.set <- as.vector(off)

      if(object$var.prop){
        var.prop <- dsm.var.prop(object$dsm.obj,object$pred.data,object$off.set,
                                 object$seglen.varname, object$type.pred)
      }else{
        var.prop <- dsm.var.gam(object$dsm.obj,object$pred.data,object$off.set,
                                 object$seglen.varname, object$type.pred)
      }

      sinfo$se <- sqrt(var.prop$pred.var)
    }
    # grab the predicted values
#    if(length(object$pred)==0){
#      dsm.obj <- object$dsm.obj
#      class(dsm.obj) <- c("dsm",class(dsm.obj))
#      mod1.pred <- predict(dsm.obj,
#                           newdata=object$pred.data,
#                           off.set=object$off.set)
#      sinfo$pred.est <- sum(mod1.pred,na.rm=TRUE)
#    }else if(length(object$pred)>1 | length(object$pred)==0){
    if(length(object$pred)>1){
      sinfo$pred.est <- sum(unlist(object$pred), na.rm=TRUE)
    }else{
      sinfo$pred.est <- object$pred[[1]]
    }

    # if we're using variance propagation or there is no detection
    # function, then the CV is fine
    if(sinfo$varprop | is.null(object$dsm.object$ddf)){
      # calculate the CV
      sinfo$cv <- sinfo$se/sinfo$pred.est
    }else{
    # if we're just using the GAM variance, then we need to combine using
    # the delta method
      ddf.summary<-summary(object$dsm.object$ddf)

      cvp.sq <- (ddf.summary$average.p.se/
                 ddf.summary$average.p)^2
      sinfo$detfct.cv <- sqrt(cvp.sq)
      sinfo$gam.cv <- sinfo$se/sinfo$pred.est

      sinfo$cv <- sqrt(cvp.sq+sinfo$gam.cv^2)
    }
  }

  class(sinfo) <- "summary.dsm.var"
  return(sinfo)
}
