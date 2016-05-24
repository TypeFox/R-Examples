
################################################################################
## effic.make.link provides the linkinv and the mu.eta functions which are also
## provided by make.link or if applicable make.vslink
## (these two functions additionally pass back some other things)
## for certain links however, mu.eta and linkinv as passed by make.link
## is computationally quite unefficient, we provide here the more
## efficient functions in these cases
## the advantage is that the family object itself does not need to be changed
################################################################################

effic.make.link <- function(link)  {

    if(link=="varstab"){
       expr <- attr(link, "expression")
       linkvarstab <- TRUE
       link <- attr(link, "name")
       if(is.null(link)) link <- ""
    } else linkvarstab <- FALSE

    if(link=="log") {   # this is more efficient than pmax(exp(eta), .Machine$double.eps)
        linkinv.fct <- mu.eta.fct <-   function(eta){
        pm <- exp(eta)
        pm[pm<.Machine$double.eps] <- .Machine$double.eps
        pm  }
    } else if(link=="probit" ){
        mu.eta.fct <- function(eta){
            pm <- dnorm(eta)
            pm[pm<.Machine$double.eps] <- .Machine$double.eps
            pm  }
        linkinv.fct <- function(eta) {
            thresh <- -qnorm(.Machine$double.eps)
		        # eta <- pmax(eta, -thresh)
		        eta[eta < - thresh] <- - thresh
            # eta <- pmin(eta , thresh)
		        eta[eta > thresh] <- thresh
            pnorm(eta)         }
    } else if(link=="cauchit" ){
         mu.eta.fct <- function(eta){
            pm <- dcauchy(eta)
              pm[pm<.Machine$double.eps] <- .Machine$double.eps
              pm             }
         linkinv.fct <- function(eta) {
            thresh <- -qcauchy(.Machine$double.eps)
		        # eta <- pmax(eta, -thresh)
		        eta[eta < - thresh] <- - thresh
            # eta <- pmin(eta , thresh)
		        eta[eta > thresh] <- thresh
            pcauchy(eta)            }
    } else if(link=="cloglog" ){
         mu.eta.fct <- function(eta){
            # eta <- pmin(eta, 700)
            eta[eta>700] <- 700
            pm <- exp(eta) * exp(-exp(eta))
            pm[pm<.Machine$double.eps] <- .Machine$double.eps
            pm               }
        linkinv.fct <- function(eta) {
            pm <- -expm1(-exp(eta))
            # pm <- pmin( pm , 1 - .Machine$double.eps)
            pm[pm >  1 - .Machine$double.eps] <-  1 - .Machine$double.eps
            # pmax(pm, .Machine$double.eps)
            pm[pm < .Machine$double.eps] <- .Machine$double.eps             }
    } else {
        if(linkvarstab & link=="" ){ # is.expression(expr) # attr(link, "name")
            # there are also varstab link with attr like "mu/(mu-1)" that work with make.link
              mk <- eval( expr)
              # mu.eta.fct <- family$mu.eta
              # linkinv.fct <- family$linkinv
        } else {

          mk <- make.link(link)
        }
          mu.eta.fct  <- mk$mu.eta
          linkinv.fct <- mk$linkinv

    }
    return(list(linkinv.fct=linkinv.fct, mu.eta.fct=mu.eta.fct) )

}

