#' Variance estimation via parametric moving block bootstrap
#'
#' Estimate the variance in abundance over an area using a moving block
#' bootstrap. Two procedures are implemented, one incorporating detection
#' function uncertainty, one not.
#'
#' @inheritParams dsm.var.gam
#' @param dsm.object object returned from \code{\link{dsm}}.
#' @param n.boot number of bootstrap resamples.
#' @param block.size number of segments in each block.
#' @param ds.uncertainty incorporate uncertainty in the detection function? See
#'        Details, below. Note that this feature is EXPERIMENTAL at the moment.
#' @param samp.unit.name name sampling unit to resample (default
#'        'Transect.Label').
#' @param progress.file path to a file to be used (usually by Distance) to
#'        generate a progress bar (default \code{NULL} -- no file written).
#' @param bs.file path to a file to store each boostrap round. This stores all
#'        of the bootstrap results rather than just the summaries, enabling
#'        outliers to be detected and removed. (Default \code{NULL}).
#' @param bar should a progress bar be printed to screen? (Default \code{TRUE}).
#'
#' @section Details:
#'  Setting \code{ds.uncertainty=TRUE} will incorporate detection function
#'  uncertainty directly into the bootstrap. This is done by generating
#'  observations from the fitted detection function and then re-fitting a new
#'  detection function (of the same form), then calculating a new effective
#'  strip width. Rejection sampling is used to generate the observations (except
#'  in the half-normal case) so the procedure can be rather slow. Note that this
#'  is currently not supported with covariates in the detection function.
#'
#'  Setting \code{ds.uncertainty=FALSE} will incorporate detection function
#'  uncertainty using the delta method. This assumes that the detection function
#'  and the spatial model are INDEPENDENT. This is probably not reasonable.
#'
#' @export
#' @importFrom utils write.table setTxtProgressBar txtProgressBar
#' @examples
#' \dontrun{
#' library(Distance)
#' library(dsm)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
#'                key = "hr", adjustment = NULL)
#' summary(hr.model)
#'
#' # fit a simple smooth of x and y
#' mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
#' summary(mod1)
#'
#' # create an offset (in metres)
#' # each prediction cell is 444km2
#' off.set <- 444*1000*1000
#'
#' # calculate the variance by 500 moving block bootstraps
#' mod1.movblk <- dsm.var.movblk(mod1, mexdolphins$preddata, n.boot = 500,
#'    block.size = 3, samp.unit.name = "Transect.Label", off.set = off.set,
#'    bar = TRUE, bs.file = "mexico-bs.csv", ds.uncertainty = TRUE)
#' }

## TODO
# make up the missed replicates due to model errors?
# detection function uncertainty:
#  * new sampler doesn't do point transects at the moment
#  * new sampler doesn't do covariates at the moment
#  * is the sample size correct?
#  * should we be calculating the offset in here?

# this used to be called param.movblk.variance
dsm.var.movblk <- function(dsm.object, pred.data, n.boot, block.size,
                           off.set, ds.uncertainty=FALSE,
                           samp.unit.name='Transect.Label',
                           progress.file=NULL, bs.file=NULL,bar=TRUE){

  # check the user didn't ask for DS uncertainty and didn't supply
  # a detection function
  if(ds.uncertainty & is.null(dsm.object$ddf)){
    stop("Cannot incorporate detection function uncertainty with no detection function!")
  }
  # check the user didn't ask for individual level covars and detection
  # function uncertainty
  if(ds.uncertainty & !is.null(dsm.object$ddf) &&
     dsm.object$ddf$ds$aux$ddfobj$scale$formula != "~1"){
    stop("Cannot incorporate detection function uncertainty with covariates in the detection function")
  }

  # Initialize storage
  study.area.total <- numeric(n.boot)
  short.var <- data.frame(sumx=rep(0,nrow(pred.data)),
                          sumx.sq=rep(0,nrow(pred.data)))

  # save the original off.set that was supplied
  original.offset<-off.set
  # extract the link and inverse link functions
  linkfn <- dsm.object$family$linkfun
  invlinkfn <- dsm.object$family$linkinv
  # what was the response?
  response <- as.character(dsm.object$formula)[2]

  # Sort out sampling unit for dsm object
  dsm.object$data$sampling.unit <- dsm.object$data[[samp.unit.name]]

  # Following line removes any transect into which missing data was
  # detected by the call to gam and recorded in '$na.action'
  # Consequence of this step should be that no sampling unit (transect)
  # can become part of the bootstrap when any segment has missing data.
  name.sampling.unit <- unique(dsm.object$data$sampling.unit)
  num.sampling.unit <- length(name.sampling.unit)

  # Get residuals
  obs <- dsm.object$data[[response]]
  fit.vals <- rep(NA,length(obs))
  fit.vals[!is.na(obs)] <- fitted(dsm.object)
  eps <- 0.001
  dsm.object$data$log.resids <- linkfn(obs+eps)-linkfn(fit.vals+eps)

  # Sort out blocks for each sampling unit
  block.info <- block.info.per.su(block.size=block.size,
                                  data=dsm.object$data,
                                  name.su=name.sampling.unit)
  tot.num.blocks <- sum(block.info$num.block)
  num.blocks.required <- sum(block.info$num.req)
  block.vector <- 1:tot.num.blocks

  # do we want to print a progress bar?
  if(bar){
    pb <- txtProgressBar(min=0,max=n.boot,style=3)
  }

  # first resample is the model predicitons themselves

  # if we supplied bs.file, then write detailed, per replicate
  # if it's the first time and the file already exists, modify the name 
  # use that
  if(!is.null(bs.file)){
    if(file.exists(bs.file)){
      a.number <- 1
      bs.file2 <- strsplit(bs.file,"\\.")[[1]]
      bs.file2.end <- paste('.',bs.file2[length(bs.file2)],sep="")
      bs.file2.start <- paste(bs.file2[1:(length(bs.file2)-1)],collapse=".")
      bs.file2 <- paste(bs.file2.start,"-",a.number,bs.file2.end,
                        collapse="",sep="")
      while(file.exists(bs.file2)){
        a.number <- a.number+1
        bs.file2 <- paste(bs.file2.start,"-",a.number,bs.file2.end,
                          collapse="",sep="")
      }
      warning(paste("Filename",bs.file,"was taken, writing to file",
                     bs.file2))
      bs.file<-bs.file2
    }
  }

  ## prediction from the observed data
  dsm.predict.bs <- try(predict(dsm.object,
                                newdata=pred.data,
                                off.set=off.set))
  # we shouldn't get errors in this predict, but check just in case
  if(all(class(dsm.predict.bs)=="try-error")){
    dsm.predict.bs <- rep(NA,length(short.var$sumx))
  }

  # Don't save all cell values for all reps, rather,
  #  populate dataframe with machine formula components for each cell
  short.var$sumx <- short.var$sumx + dsm.predict.bs
  short.var$sumx.sq <- short.var$sumx.sq + dsm.predict.bs^2
  study.area.total[1] <- sum(dsm.predict.bs, na.rm=TRUE)

  if(!is.null(bs.file)){
    # append to the file bs.file
    write.table(t(dsm.predict.bs),bs.file,append=TRUE,
                sep=",",col.names=FALSE)
  }

  # save the data used to fit the model
  data.save <- dsm.object$data

  # Start bootstrapping
  for(i in 2:n.boot){
    # Compute proportion of bootstrapping completed, write it to file
    # passed as argument 'progress.file'
    # This file will be read by Distance to present a progress bar.
    if(!is.null(progress.file)){
      progress <- round(i/n.boot, 2)* 100
      write(progress, file=progress.file, append=FALSE)
    }

    # make the boostrap sample of the blocks
    bs.blocks <- sample(block.vector, num.blocks.required, replace=TRUE)
    # extract the sample residuals from the data
    bs.resids <- generate.mb.sample(num.blocks.required, block.size,
                                    bs.blocks, dsm.object$data,
                                    block.info, num.sampling.unit)
    # copy the data into a new object
    bs.samp <- data.save

    # if we are incorporating detection function uncertainty, then
    # need to resample the distances
    if(ds.uncertainty){

      # save the old probability of detection
      old.p<-fitted(dsm.object$ddf)[1]

      #### This only deals with count data at the moment
      ####  => no individual level covariates
      #### Much of this can be put up top, doesn't need to be calculated
      #### each time

      # call out to generate some ds data, and fit a model to that data,
      # asumming that the detection function model is correct
      new.p <- generate.ds.uncertainty(dsm.object$ddf)

      # calculate the new offset
      this.offset <- new.p*(invlinkfn(dsm.object$data$off.set)/old.p)

      # calculate the fitted values with the new offset
      fit <- (fitted(dsm.object)/invlinkfn(dsm.object$data$off.set))*
                                    this.offset

      # replace the offset in the model
      bs.samp$off.set <- linkfn(this.offset)
    }else{
      # if we're not doing detection function uncertainty then
      #  just get the fitted values
      fit <- fitted(dsm.object)
    }

    if(any(is.na(fit))){
      stop("Data with NA covariate values can't be used with moving block")
    }

    # calculate the bootstrap "observations"
    bs.samp[[response]] <- fit*invlinkfn(bs.resids)


    ## Fit model to dsm bootstrap sample
    # Reconstruct dsm model fitting command -- this is a call to gam() or gamm()
    gam.call <- dsm.object$call
    gam.call$formula <- dsm.object$formula
    gam.call$family <- dsm.object$family
    ## put the bootstrap data into the gam call
    gam.call$data <- bs.samp
    # fit model to bootstrap data
    #   try() to handle model failure,
    #   withCallingHandlers() to deal with possible spurious error message
    #   with() is not required for gam but IS necessary for gamm models
    #     to obtain covariance structure etc
    #dsm.bs.model <- try(with(dsm.object,eval(gam.call)))
    dsm.bs.model <- try(with(dsm.object,withCallingHandlers(eval(gam.call),
                                           warning = matrixnotposdef.handler)))

    # if everything didn't fail...
    if(all(class(dsm.bs.model)!="try-error")){
      # make the model a dsm
      class(dsm.bs.model) <- c("dsm",class(dsm.bs.model))

      # Do prediction using newly fitted dsm model created from bootstrap sample
      dsm.predict.bs <- try(predict(dsm.bs.model,
                                    newdata=pred.data,
                                    off.set=off.set))
      if(all(class(dsm.predict.bs)=="try-error")){
        dsm.predict.bs <- rep(NA,length(fit))
      }

      # Don't save all cell values for all reps, rather,
      #  populate data.frame with formula components for each cell
      dsm.predict.bs <- dsm.predict.bs
      short.var$sumx <- short.var$sumx + dsm.predict.bs
      short.var$sumx.sq <- short.var$sumx.sq + dsm.predict.bs^2
      study.area.total[i] <- sum(dsm.predict.bs, na.rm=TRUE)

      if(!is.null(bs.file)){
        # append predictions to the file bs.file
        write.table(t(dsm.predict.bs),bs.file,append=TRUE,
                    sep=",",col.names=FALSE)
      }
    }else{
      study.area.total[i] <- NA
    }

    # the model and data might be quite big, so clear some space
    rm(dsm.predict.bs,gam.call,bs.samp)
    gc()

    # add a tick to the progress bar
    if(bar){
      setTxtProgressBar(pb, i)
    }
  }

  result <- list(short.var=short.var,
                 study.area.total=study.area.total,
                 ds.uncertainty=ds.uncertainty,
                 bootstrap=TRUE,
                 pred.data=pred.data,
                 n.boot=n.boot,
                 off.set=original.offset,
                 block.size=block.size,
                 bs.file=bs.file,
                 dsm.object=dsm.object
                )

  # give it a bit of class
  class(result)<-c("dsm.var")

  if(bar){
    cat("\n")
  }
  return(result)
}
