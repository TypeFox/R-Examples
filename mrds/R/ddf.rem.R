#' Mark-Recapture Distance Sampling (MRDS) Removal - PI
#'
#' Mark-Recapture Distance Sampling (MRDS) Analysis of Removal Observer
#' Configuration and Point Independence
#'
#' MRDS analysis based on point independence involves two separate and
#' independent analyses of the mark-recapture data and the distance sampling
#' data.  For the removal observer configuration, the mark-recapture data are
#' analysed with a call to \code{\link{ddf.rem.fi}} (see Laake and Borchers
#' 2004) to fit conditional distance sampling detection functions to estimate
#' p(0), detection probability at distance zero for the primary observer based
#' on independence at zero (eq 6.22 in Laake and Borchers 2004). Independently,
#' the distance data, the observations from the primary observer, are used to
#' fit a conventional distance sampling (CDS) (likelihood eq 6.6) or
#' multi-covariate distance sampling (MCDS) (likelihood eq 6.14) model for the
#' detection function, g(y), such that g(0)=1. The detection function for the
#' primary observer is then created as p(y)=p(0)*g(y) (eq 6.28 of Laake and
#' Borchers 2004) from which predictions are made. \code{ddf.rem} is not called
#' directly by the user and is called from \code{\link{ddf}} with
#' \code{method="rem"}.
#'
#' For a complete description of each of the calling arguments, see
#' \code{\link{ddf}}.  The argument \code{data} is the dataframe specified by
#' the argument \code{data} in \code{ddf}. The arguments \code{dsmodel},
#' \code{mrmodel}, \code{control} and \code{meta.data} are defined the same as
#' in \code{ddf}.
#'
#' @export
#' @method ddf rem
#' @param dsmodel distance sampling model specification; model list with key
#'   function and scale formula if any
#' @param mrmodel mark-recapture model specfication; model list with formula
#'   and link
#' @param data analysis dataframe
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call original function call used to call \code{ddf}
#' @return result: an rem model object which is composed of rem.fi and ds model
#'   objects
#' @author Jeff Laake
#' @seealso \code{\link{ddf.rem.fi}}, \code{\link{ddf.ds}}
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#' @keywords Statistical Models
ddf.rem<-function(dsmodel,mrmodel,data,meta.data=list(),control=list(),call=""){


  # Test to make sure that observer not used in mrmodel
  if(length(grep("observer",mrmodel))!=0){
    stop("observer cannot be included in models for removal configurations\n")
  }
  # Save current user options and then set design contrasts to treatment style
  save.options<-options()
  options(contrasts=c("contr.treatment","contr.poly"))

  # Set up meta data values
  meta.data <- assign.default.values(meta.data, left=0, width=NA, binned=FALSE,
                                     int.range=NA, point=FALSE)

  # Set up control values
  control <- assign.default.values(control, showit=0,
                                   estimate=TRUE, refit=TRUE, nrefits=25,
                                   initial=NA, lowerbounds=NA, upperbounds=NA,
                                   mono.points=20)

  # Process data
  data.list <- process.data(data,meta.data)
  meta.data <- data.list$meta.data
  xmat <- data.list$xmat

  # Create result list
  result <- list(call=call, data=data, mrmodel=mrmodel, dsmodel=dsmodel,
                 meta.data=meta.data, control=control, method="rem")
  class(result) <- c("rem","ddf")

  #  Fit the conditional detection functions using ddf.rem.fi
  result$mr <- ddf.rem.fi(model=mrmodel, data, meta.data, control, call,
                          method="rem")

  #  Fit the unconditional detection functions using ddf.ds
  #  5/24/05 - jll add call to process.data for unique.data because it
  #  didn't handle truncation correctly - error reported by Sharon Hedley
  unique.data <- data[data$observer==1&data$detected==1,]
  obs2 <- data[data$observer==2,]
  obs1 <- data[data$observer==1,]
  #missed <- data$detected[data$observer==1]==0
  unique.data <- rbind(unique.data,obs2[obs1$detected==0,])
  unique.data$observer <- 1
  unique.data  <-  process.data(unique.data, meta.data,check=FALSE)$xmat
  result$ds <- ddf.ds(model=dsmodel,unique.data,meta.data,control,call)

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
  result$criterion<- -2*result$lnl + 2*npar

  # Get fitted values and predict abundance and its variance in covered region
  result$fitted <- predict(result,newdata=xmat,compute=TRUE)$fitted
  result$Nhat <- NCovered(result)

  # Restore user options
  options(save.options)

  # Return result
  return(result)
}
