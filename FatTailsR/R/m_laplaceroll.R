

#' @include l_estimation.R





#' @title Laplace-Gauss Normal Distribution Object
#'
#' @description
#' An object designed after regkienerLX to summarize the information related to 
#' a given dataset when the Laplace-Gauss normal distribution is applied on it.
#' 
#' @param    X       vector of quantiles.
#' @details      This function is designed after regkienerLX to provide a 
#'               similar framework.

#' @return  
#' A list with the following data.frame:
#' \itemize{
#'   \item{dfrXPn}{data.frame. X = initial quantiles. Pn = estimated normal probabilites.}
#'   \item{dfrXLn}{data.frame. X = initial quantiles. Ln = logit of estimated normal probabilites.}
#'   \item{dfrXDn}{data.frame. X = initial quantiles. Dn = estimated normal density.}
#'   \item{coefn}{numeric. The mean and the standard deviation of the dataset.}
#'   \item{dfrQnPn}{data.frame. Qn = estimated quantiles of interest. Pn = probability.} 
#'   \item{dfrQnPn}{data.frame. Qn = estimated quantiles of interest. Pn = logit of probability.}
#' }
#'
#'
#' @examples     
#' prices2returns <- function(x) { 100*diff(log(x)) }
#' CAC  <- prices2returns(as.numeric(EuStockMarkets[,3])) 
#' lgn  <- laplacegaussnorm( CAC )
#' attributes(lgn)
#' head(lgn$dfrXPn)
#' head(lgn$dfrXLn)
#' head(lgn$dfrXDn)
#' lgn$coefn
#' lgn$dfrQnPn
#' lgn$dfrQnLn
#'
#' @seealso      The regression function \code{\link{regkienerLX}}.
#' @export
#' @name laplacegaussnorm
         laplacegaussnorm <- function(X) {

X       <- sort(as.numeric(X)) 
mX      <- mean(X)
sX      <- sd(X)
proban  <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.5,
               0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
quantn  <- qnorm(p = proban, mean = mX, sd = sX )
names(quantn)  <- c("n.0001", "n.0005", "n.001", "n.005", "n.01", "n.05", 
            "n.50", "n.95", "n.99", "n.995", "n.999", "n.9995", "n.9999")
dfrQnPn <- data.frame( Qn = quantn, Pn = proban )
dfrQnLn <- data.frame( Qn = quantn, Ln = logit(proban) )

# Final objet lgn
lgn     <- list()
lgn$dfrXPn    <- data.frame( X = X, Pn = pnorm(X, mX, sX) )
lgn$dfrXLn    <- data.frame( X = X, Ln = logit(pnorm(X, mX, sX)) )
lgn$dfrXDn    <- data.frame( X = X, Dn = dnorm(X, mX, sX) )
lgn$coefn     <- c( m = mX, sd = sX )
lgn$dfrQnPn   <- dfrQnPn
lgn$dfrQnLn   <- dfrQnLn
class(lgn)    <- "cllgn"

return(lgn)
}


