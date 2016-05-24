



#' @title Package FatTailsR 
#' 
#' @description
#' This package includes Kiener distributions K1, K2, K3 and K4 and two estimation 
#' algorithms to characterize with a high precision symmetric or asymmetric 
#' distributions with left and right fat tails that appear in market finance, 
#' neuroscience and many other disciplines. The estimation of the distribution parameters, 
#' quantiles, value-at-risk and expected shortfall is usually very accurate. 
#' Two datasets are provided, as well as power hyperbolas and power hyperbolic 
#' functions which are simplified versions of symmetric distribution K1.
#' Some functions introduced in v1.2-0 were discarded or renamed.
#' 
#' Download the pdf cited in the references to get an overview of  
#' the theoretical part and several examples on stocks and indices. 
#' 
#' A commercial package, \code{FatTailsRplot}, with advanced plotting functions 
#' and calculation of matrix of stocks over rolling windows is also developped 
#' by the author. 
#' 
#' @details
#' With so many functions, this package could look fat. But it's not! 
#' It's rather agile and easy to use! The various functions included in this package 
#' can be assigned to the following groups:
#' \enumerate{
#'   \item Several vectors of probabilities: 
#'         \itemize{
#'         \item \code{\link{pprobs0}}, pprobs1, pprobs2, ..., pprobs9.
#'         }
#'   \item Two datasets presented in different formats:   
#'         list, data.frame, timeSeries, xts, zoo:
#'         \itemize{
#'         \item \code{\link{getDSdata}}.
#'         \item \code{\link{extractData}}, dfData, tData, xData, zData.
#'         }
#'   \item  Miscellaneous functions and functions related to the logistic function:
#'         \itemize{
#'         \item \code{\link{dimdim}}, dimdim1.
#'         \item \code{\link{logit}}, invlogit, ltmlogis, rtmlogis, eslogis.
#'         }
#'   \item Power hyperbolas, power hyperbolic functions and their inverses: 
#'         \itemize{
#'         \item \code{\link{exphp}}, coshp, sinhp, tanhp, sechp, cosechp, 
#'               cotanhp.
#'         \item \code{\link{loghp}}, acoshp, asinhp, atanhp, asechp, 
#'               acosechp, acotanhp.
#'         \item \code{\link{kashp}}, dkashp_dx, ashp.
#'         }
#'   \item Logishp function, kogit and invkogit = logistic function + power hyperbolas: 
#'         \itemize{
#'         \item d, p, q, r, dp, dq, l, dl, ql \code{\link{logishp}}.
#'         \item \code{\link{kogit}}, invkogit.
#'         }
#'   \item Conversion functions between parameters of Kiener distributions K1, K2, K3, K4:
#'         \itemize{
#'         \item \code{\link{aw2k}}, aw2d, aw2e, ad2e, ad2k, ad2w, ae2d, ae2k, 
#'                 ae2w, ak2e, ak2w, de2a, de2k, de2w, dk2a, dk2e, dw2a, dw2e, 
#'                 dw2k, ek2a, ak2d, ek2w, aw2a, aw2d, ew2a, aw2d, ew2k, kd2a, 
#'                 kd2e, kd2w, ke2a, ke2d, ke2w, kw2a, kw2d, kw2e.
#'         \item \code{\link{pk2pk}}.
#'         }
#'   \item Kiener distributions K1, K2, K3, K4:
#'         \itemize{
#'         \item d, p, q, r, dp, dq, l, dl, ql, var, ltm, rtm, dtmq, es \code{\link{kiener1}},
#'         \item d, p, q, r, dp, dq, l, dl, ql, var, ltm, rtm, dtmq, es \code{\link{kiener2}},
#'         \item d, p, q, r, dp, dq, l, dl, ql, var, ltm, rtm, dtmq, es \code{\link{kiener3}},
#'         \item d, p, q, r, dp, dq, l, dl, ql, var, ltm, rtm, dtmq, es \code{\link{kiener4}}.
#'         }
#'   \item Quantile (VaR) corrective function (as a multiplier of the logistic function). 
#'         Expected shortfall corrective function (as a multiplier of the expected shortfall 
#'         of the logistic distribution):
#'         \itemize{
#'         \item \code{\link{ckiener1}}, ckiener2, ckiener3, ckiener4.
#'         \item \code{\link{hkiener1}}, hkiener2, hkiener3, hkiener4.
#'         }
#'   \item Moments of the distribution estimated from the dataset and from 
#'         the regression parameters:
#'         \itemize{
#'         \item \code{\link{xmoments}}.
#'         \item \code{\link{kmoments}}, kmoment, kcmoment, kmean, 
#'               kstandev, kvariance, kskewness, kkurtosis, kekurtosis.
#'         }
#'   \item Regression and estimation functions to estimate Kiener distribution
#'         parameters on a given dataset. \code{*fit*} and \code{*param*} 
#'         are wrappers of algorithms \code{reg} and \code{estim}.
#'         \code{reg} uses an unweighted nonlinear regression function. 
#'         \code{estim} uses a fast estimation based on quantiles:
#'         \itemize{
#'         \item \code{\link{regkienerLX}}, \code{\link{laplacegaussnorm}}.
#'         \item \code{\link{fitkienerX}}.
#'         \item \code{\link{paramkienerX}, paramkienerX5, paramkienerX7}.
#'         }
#'   \item Functions related to \code{paramkienerX}:
#'         \itemize{
#'         \item \code{\link{elevenprobs}}, sevenprobs, fiveprobs.
#'         \item \code{\link{estimkiener11}}, estimkiener7, estimkiener5.
#'         \item \code{\link{roundcoefk}}.
#'         }
#'   \item Predefined parameter subsets to extract the corresponding parameters 
#'         from the long vector \code{fitk} obtained after regression/estimation 
#'         \code{regkienerLX}, \code{fitkienerX} :
#'         \itemize{
#'         \item \code{\link{exfit0}}, ..., \code{exfit7}.
#'         }
#' }
#' For a quick start, jump to the functions \code{\link{regkienerLX}},  
#' \code{\link{fitkienerX}} and run the examples. 
#' Then, download and read the documents in pdf format cited in the references 
#' to get an overview on the major functions. Finally, explore the other 
#' examples. 
#' 
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014. Download it from: 
#' \url{http://www.inmodelia.com/exemples/2014-0627-Rmetrics-Kiener-en.pdf}
#'
#' P. Kiener, Fat tail analysis and package FatTailsR, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' Download it from: 
#' \url{http://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#'
#' @keywords symbolmath distribution models
#' 
#' @examples     
#' 
#' require(graphics)
#' require(minpack.lm)
#' require(timeSeries)
#' 
#' ### Load the datasets and select one number (1-16)
#' DS     <- getDSdata()
#' j      <- 5
#' 
#' ### and run this block
#' X      <- DS[[j]]
#' nameX  <- names(DS)[j]
#' reg    <- regkienerLX(X)
#' lgn    <- laplacegaussnorm(X)
#' lleg   <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", 
#'            "logit(0.95)   = 2.9", "logit(0.50)   = 0", 
#'            "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'            "logit(0.001) = -6.9  ")
#' pleg   <- c( paste("m =",  reg$coefr4[1]), paste("g  =", reg$coefr4[2]), 
#'              paste("k  =", reg$coefr4[3]), paste("e  =", reg$coefr4[4]) )
#' 
#' ## Main plot
#' op     <- par(mfrow = c(1,1), mgp = c(1.5,0.8,0), mar = c(3,3,2,1))
#' plot(reg$dfrXP, main = nameX)
#' legend("top", legend = pleg, cex = 0.9, inset = 0.02 )
#' lines(reg$dfrEP, col = 2, lwd = 2)
#' points(reg$dfrQkPk, pch = 3, col = 2, lwd = 2, cex = 1.5)
#' lines(lgn$dfrXPn, col = 7, lwd = 2)
#' 
#' ## Plot F(X) > 0,97
#' front = c(0.06, 0.39, 0.50, 0.95)
#' par(fig = front, new = TRUE, mgp = c(1.5, 0.6, 0), las = 0)
#' plot( reg$dfrXP[which(reg$dfrXP$P > 0.97),] , pch = 1, xlab = "", ylab = "", main = "F(X) > 0,97" )
#' lines(reg$dfrEP[which(reg$dfrEP$P > 0.97),], type="l", col = 2, lwd = 3 )
#' lines(lgn$dfrXPn[which(lgn$dfrXPn$Pn > 0.97),], type = "l", col = 7, lwd= 2 )
#' points(reg$dfrQkPk, pch = 3, col = 2, lwd = 2, cex = 1.5)
#' points(lgn$dfrQnPn, pch = 3, col = 7, lwd = 2, cex = 1)
#' 
#' ## Plot F(X) < 0,03
#' front = c(0.58, 0.98, 0.06, 0.61)
#' par(fig = front, new = TRUE, mgp = c(0.5, 0.6, 0), las = 0 )
#' plot( reg$dfrXP[which(reg$dfrXP$P < 0.03),] , pch = 1, xlab = "", ylab = "", main = "F(X) < 0,03")
#' lines(reg$dfrEP[which(reg$dfrEP$P < 0.03),], type = "l", col = 2, lwd = 3 )
#' lines(lgn$dfrXPn[which(lgn$dfrXPn$Pn < 0.03),], type = "l", col= 7, lwd= 2 )
#' points(reg$dfrQkPk, pch = 3, col = 2, lwd = 2, cex = 1.5)
#' points(lgn$dfrQnPn, pch = 3, col = 7, lwd = 2, cex = 1)
#' 
#' ## Moments from the parameters (k) and from the Dataset (X)
#' round(cbind("k" = kmoments(reg$coefk, lengthx = nrow(reg$dfrXL)), "X" = xmoments(X)), 2)
#' attributes(reg)
#' ### End block
#' 
#' 
#' @importFrom  parallel detectCores makeCluster stopCluster parApply mclapply
#' @import  timeSeries
#' @import  minpack.lm
#' @import  stats
#' @rdname  FatTailsR
#' @aliases FatTailsR
#' @name    FatTailsR-package
#' @docType package
NULL


