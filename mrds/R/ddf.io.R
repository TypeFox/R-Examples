#' Mark-Recapture Distance Sampling (MRDS) IO - PI
#'
#' Mark-Recapture Distance Sampling (MRDS) Analysis of Independent Observer
#' Configuration and Point Independence
#'
#' MRDS analysis based on point independence involves two separate and
#' independent analyses of the mark-recapture data and the distance sampling
#' data.  For the independent observer configuration, the mark-recapture data
#' are analysed with a call to \code{\link{ddf.io.fi}} (see likelihood eq 6.8
#' and 6.16 in Laake and Borchers 2004) to fit conditional distance sampling
#' detection functions to estimate p(0), detection probability at distance zero
#' for the independent observer team based on independence at zero (eq 6.22 in
#' Laake and Borchers 2004). Independently, the distance data, the union of the
#' observations from the independent observers, are used to fit a conventional
#' distance sampling (CDS) (likelihood eq 6.6) or multi-covariate distance
#' sampling (MCDS) (likelihood eq 6.14) model for the detection function, g(y),
#' such that g(0)=1. The detection function for the observer team is then
#' created as p(y)=p(0)*g(y) (eq 6.28 of Laake and Borchers 2004) from which
#' predictions are made. \code{ddf.io} is not called directly by the user and
#' is called from \code{\link{ddf}} with \code{method="io"}.
#'
#' For a complete description of each of the calling arguments, see
#' \code{\link{ddf}}.  The argument \code{dataname} is the name of the
#' dataframe specified by the argument \code{data} in \code{ddf}. The arguments
#' \code{dsmodel}, \code{mrmodel}, \code{control} and \code{meta.data} are
#' defined the same as in \code{ddf}.
#'
#' @export
#' @method ddf io
#' @param dsmodel distance sampling model specification; model list with key
#'   function and scale formula if any
#' @param mrmodel mark-recapture model specfication; model list with formula
#'   and link
#' @param data analysis dataframe
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call original function call used to call \code{ddf}
#' @return result: an io model object which is composed of io.fi and ds model
#'   objects
#' @author Jeff Laake
#' @seealso \code{\link{ddf.io.fi}},
#'   \code{\link{ddf.ds}},\code{\link{summary.io}},\code{\link{coef.io}},\code{\link{plot.io}},
#'   \code{\link{gof.io}}
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#' @keywords Statistical Models
ddf.io<-function(dsmodel,mrmodel,data,meta.data=list(),control=list(),call=""){

  # Save current user options and then set design contrasts to treatment style
  save.options<-options()
  options(contrasts=c("contr.treatment","contr.poly"))

  # Set up meta data values
  meta.data=assign.default.values(meta.data, left=0, width=NA, binned=FALSE,
                                   int.range=NA,point=FALSE)

  # Set up control values
  control=assign.default.values(control, showit=0,
                                estimate=TRUE, refit=TRUE, nrefits=25,
                                initial=NA, lowerbounds=NA, upperbounds=NA,
                                mono.points=20)

  # Process data
  data.list <- process.data(data,meta.data)
  meta.data <- data.list$meta.data

  # Create result list
  result <- list(call=call, data=data, mrmodel=mrmodel, dsmodel=dsmodel,
                 meta.data=meta.data, control=control, method="io")
  class(result)=c("io","ddf")

  # Fit the conditional detection functions using ddf.io.fi  
  result$mr <- ddf.io.fi(model=mrmodel,data,meta.data,control,call,method="io")

  # Fit the unconditional detection functions using ddf.ds
  unique.data <- data[data$observer==1,]
  unique.data$detected <- 1
  result$ds <- ddf.ds(model=dsmodel,unique.data,meta.data,control,call)
  # use the ds meta data (as this has correct mono values)
  result$meta.data <- result$ds$meta.data
  if(is.null(result$ds$Nhat)){
    if(control$debug){
      warning("ds model did not converge; no further results possible\nReturned object is for debugging ONLY!")
      return(result)
    }
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
  result$fitted=predict(result)$fitted
  result$Nhat=NCovered(result)

  # Restore user options
  options(save.options)

  # Return result
  return(result)
}
