#' Mark-Recapture Distance Sampling (MRDS) Trial Configuration - PI
#'
#' Mark-Recapture Distance Sampling (MRDS) Analysis of Trial Observer
#' Configuration and Point Independence
#'
#' MRDS analysis based on point independence involves two separate and
#' independent analyses of the mark-recapture data and the distance sampling
#' data.  For the trial configuration, the mark-recapture data are analysed
#' with a call to \code{\link{ddf.trial.fi}} (see likelihood eq 6.12 and 6.17
#' in Laake and Borchers 2004) to fit a conditional distance sampling detection
#' function for observer 1 based on trials (observations) from observer 2 to
#' estimate p_1(0), detection probability at distance zero for observer 1.
#' Independently, the distance data from observer 1 are used to fit a
#' conventional distance sampling (CDS) (likelihood eq 6.6) or multi-covariate
#' distance sampling (MCDS) (likelihood eq 6.14) model for the detection
#' function, g(y), such that g(0)=1.  The detection function for observer 1 is
#' then created as p_1(y)=p_1(0)*g(y) (eq 6.28 of Laake and Borchers 2004) from
#' which predictions are made. \code{ddf.trial} is not called directly by the
#' user and is called from \code{\link{ddf}} with \code{method="trial"}.
#'
#' For a complete description of each of the calling arguments, see
#' \code{\link{ddf}}.  The argument \code{dataname} is the name of the
#' dataframe specified by the argument \code{data} in \code{ddf}. The arguments
#' \code{dsmodel}, \code{mrmodel}, \code{control} and \code{meta.data} are
#' defined the same as in \code{ddf}.
#'
#' @export
#' @method ddf trial
#' @param dsmodel distance sampling model specification; model list with key
#'   function and scale formula if any
#' @param mrmodel mark-recapture model specfication; model list with formula
#'   and link
#' @param data analysis dataframe
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call original function call used to call \code{ddf}
#' @return result: a trial model object which is composed of trial.fi and ds
#'   model objects
#' @author Jeff Laake
#' @seealso \code{\link{ddf.trial.fi}},
#'   \code{\link{ddf.ds}},\code{\link{summary.trial}},\code{\link{coef.trial}},\code{\link{plot.trial}},
#'   \code{\link{gof.trial}}
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#' @keywords Statistical Models
ddf.trial <- function(dsmodel,mrmodel,data,meta.data=list(),control=list(),
                      call=""){


  # Test to make sure that observer not used in mrmodel
  if(length(grep("observer",mrmodel))!=0){
    stop("observer cannot be included in models for trial configurations\n")
  }

  # Save current user options and then set design contrasts to treatment style
  save.options <- options()
  options(contrasts=c("contr.treatment","contr.poly"))

  # Set up meta data values
  meta.data <- assign.default.values(meta.data, left=0, width=NA, binned=FALSE,
                                     int.range=NA, point=FALSE)

  # Set up control values
  control <- assign.default.values(control, showit = 0,
                                   estimate=TRUE, refit=TRUE, nrefits=25,
                                   initial=NA, lowerbounds=NA, upperbounds=NA,
                                   mono.points=20)

  # Process data
  data.list <- process.data(data,meta.data)
  meta.data <- data.list$meta.data
  xmat <- data.list$xmat

  # Create result list
  result=list(call=call, data=data, mrmodel=mrmodel, dsmodel=dsmodel,
              meta.data=meta.data, control=control, method="trial")
  class(result) <- c("trial","ddf")

  # Fit the conditional detection functions using ddf.trial.fi
  result$mr <- ddf.trial.fi(model=mrmodel,data,meta.data,control,
                            call,method="trial")

  #  Fit the unconditional detection functions using ddf.ds
  #  5/24/05 - jll add call to process.data for unique.data because it
  #  didn't handle truncation correctly - error reported by Sharon Hedley
  unique.data <- data[data$observer==1 & data$detected==1,]
  unique.data <- process.data(unique.data, meta.data,check=FALSE)$xmat
  result$ds <- ddf.ds(model=dsmodel, unique.data, meta.data, control, call)

  # stop if ds model didn't converge
  if(is.null(result$ds$Nhat)){
    stop("ds model did not converge; no further results possible")
  }

  # Combine parameter vectors and hessian matrices 
  npar.uncond <- length(result$ds$par)
  npar <- npar.uncond+length(result$mr$par)
  hessian1 <- result$mr$hessian
  if(npar.uncond==0){
    result$hessian <- hessian1
  }else{
    hessian1 <- cbind(hessian1,matrix(0,ncol=npar.uncond,nrow=npar-npar.uncond))
    hessian2 <- cbind(matrix(0,ncol=npar-npar.uncond,nrow=npar.uncond),
                      result$ds$hessian)
    result$hessian <- rbind(hessian1,hessian2)
  }

  result$par <- coef(result)
  row.names(result$hessian) <- row.names(result$par)
  colnames(result$hessian) <- row.names(result$par)
  result$par <- result$par$estimate
  names(result$par) <- row.names(result$hessian)

  # Compute total likelihood and AIC
  result$lnl <- result$ds$lnl + result$mr$lnl
  result$criterion <- -2*result$lnl + 2*npar

  # Get fitted values and predict abundance and its variance in covered region
  result$fitted <- predict(result,newdata=unique.data,compute=TRUE)$fitted
  result$Nhat <- NCovered(result)

  # Restore user options
  options(save.options)

  # Return result
  return(result)
}
