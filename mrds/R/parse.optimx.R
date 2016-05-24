#' Parse optimx results and present a nice object
#'
#' Take the resulting object from a call to optimx and make it into an object that mrds wants to talk to.
#' @param lt an optimx object
#' @return \code{lt} object that can be used later on
#' @param lnl.last last value of the log likelihood
#' @param par.last last value of the parameters
parse.optimx <- function(lt, lnl.last, par.last){
  if(any(class(lt)=="try-error") || any(is.na(lt[,1:attr(lt,"npar")]))){
    lt <- list()
    lt$conv <- 9
    lt$value <- lnl.last
    lt$par <- par.last

  }else{
    topfit.par <- coef(lt, order="value")[1, ]
    details <- attr(lt,"details")[1,]
    lt <- as.list(summary(lt, order="value")[1, ])
    lt$par <- topfit.par
    lt$message <- ""
    names(lt)[names(lt)=="convcode"] <- "conv"
    lt$hessian <- details$nhatend
  }
  return(lt)
}
