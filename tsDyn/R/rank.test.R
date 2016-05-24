#############################################
############ Rank test
#############################################



#'Test of the cointegrating rank
#'
#'Maximum-likelihood test of the cointegrating rank.
#'
#'This function computes the two maximum-likelihood tests for the cointegration
#'rank from Johansen (1996). Tests are: \describe{ \item{trace}{Test the
#'hypothesis of rank \sQuote{h} against rank \sQuote{K}, i.e. against the
#'alternative that the system is stationary.} \item{eigenvalue}{Test the
#'hypothesis of rank \sQuote{h} against rank \sQuote{h+1}.} }
#'
#'The test works for five specifications of the deterministic terms as in
#'Doornik et al (1998), to be specified in the previous call to
#'\code{\link{VECM}}: \describe{ \item{H_ql}{Unrestricted constant and trend:
#'use \code{include="both"} } \item{H_l}{Unrestricted constant and restricted
#'trend: use \code{include="const"}} and \code{LRinclude="trend"}
#'\item{H_lc}{Unrestricted constant and no trend: use \code{include="const"}}
#'\item{H_c}{Restricted constant and no trend: use \code{LRinclude="const"}}
#'\item{H_z}{No constant nor trend: use \code{include="none"}} }
#'
#'Two testing procedures can be used: \describe{ \item{Specific test}{By
#'specifying a value for \sQuote{r_null}. The \sQuote{pval} value returned
#'gives the speciifc p-value.} \item{Automatic test}{If not value is specified
#'for \sQuote{r_null}, the function makes a simple automatic test: returns the
#'rank (slot \sQuote{r}) of the first test not rejected (level specified by arg
#'\sQuote{cval}) as recommend i.a. in Doornik et al (1998, p. 544).} }
#'
#'A full table with both test statistics ad their respective p-values is given
#'in the summary method.
#'
#'P-values are obtained from the gamma aproximation from Doornik (1998, 1999).
#'Small sample values adjusted for the sample site are also available in the
#'summary method.  Note that the \sQuote{effective sample size} for the these
#'values is different from output in gretl for example.
#'
#'@aliases rank.test print.rank.test summary.rank.test
#'@param vecm \sQuote{VECM} object computed with the function
#'\code{\link{VECM}}.
#'@param type Type of test, either 'trace' or 'eigenvalue'. See details below.
#'@param r_null Rank to test specifically.
#'@param cval Critical value level for the automatic test.
#'@param x The output from \code{rank.test} for the print method.
#'@param object The output from \code{rank.test} for the summary method.
#'@param ... Unused.
#'@return An object of class \sQuote{rank.test}, with \sQuote{print} and
#'\sQuote{summary methods}.
#'@section Comparison with urca: While \code{\link[urca]{ca.jo}} in package
#'\pkg{urca} and \code{rank.test} both implement Johansen tests, there are a
#'few differences:
#'
#'\itemize{ \item \code{rank.test} gives p-values, while \code{ca.jo} gives
#'only critical values.  \item \code{rank.test} allows for five different
#'specifications of deterministic terms (see above), \code{ca.jo} for only
#'three.  \item \code{ca.jo} allows for seasonal and exogenous regressors,
#'which is not available in \code{rank.test}.  \item The lag is specified
#'differently: \code{K} from \code{ca.jo} corresponds to \code{lag}+1 in
#'\code{rank.test}.  }
#'@author Matthieu Stigler
#'@export
#'@seealso \code{\link{VECM}} for estimating a VECM. \code{\link{rank.select}}
#'to estimate the rank based on information criteria.
#'
#'\code{\link[urca]{ca.jo}} in package \pkg{urca} for another implementation of
#'Johansen cointegration test (see section \sQuote{Comparison with urca} for
#'more infos).
#'@references - Doornik, J. A. (1998) Approximations to the Asymptotic
#'Distributions of Cointegration Tests, Journal of Economic Surveys, 12, 573-93
#'
#'- Doornik, J. A. (1999) Erratum [Approximations to the Asymptotic
#'Distribution of Cointegration Tests], Journal of Economic Surveys, 13, i
#'
#'- Doornik, Hendry and Nielsen (1998) Inference in Cointegrating Models: UK M1
#'Revisited, Journal of Economic Surveys, 12, 533-72
#'
#'- Johansen, S. (1996) Likelihood-based inference in cointegrated Vector
#'Autoregresive Models, Oxford University Press
#'@keywords ts
#'@examples
#'
#'
#'data(barry)
#'
#'## estimate the VECM with Johansen! 
#'ve <- VECM(barry, lag=1, estim="ML")
#'
#'## specific test:
#'ve_test_spec <- rank.test(ve, r_null=1)
#'ve_test_spec_tr <- rank.test(ve, r_null=1, type="trace")
#'
#'ve_test_spec
#'ve_test_spec_tr
#'
#'## No specific test: automatic method
#'ve_test_unspec <- rank.test(ve)
#'ve_test_unspec_tr <- rank.test(ve, type="trace")
#'
#'ve_test_unspec
#'ve_test_unspec_tr
#'
#'## summary method: output will be same for all types/test procedure:
#'summary(ve_test_unspec_tr)
#'
#'## The function works for many specification of the VECM(), try:
#'rank.test(VECM(barry, lag=3, estim="ML"))
#'rank.test(VECM(barry, lag=3, include="both",estim="ML"))
#'rank.test(VECM(barry, lag=3, LRinclude="const",estim="ML"))
#'
#'## Note that the tests are simple likelihood ratio, and hence can be obtained also manually:
#'-2*(logLik(ve, r=1)-logLik(ve, r=2)) # eigen test, 1 against 2
#'-2*(logLik(ve, r=1)-logLik(ve, r=3)) # eigen test, 1 against 3
#'
rank.test <- function(vecm, type=c("eigen","trace"), r_null, cval=0.05){

  type <- match.arg(type)
  if(vecm$model.specific$estim!="ML") stop("Please note the 'vecm' object should be estimated with estim=' ML', not (default) OLS")

  t <- vecm$t
  k <- vecm$k
  inc <- vecm$include
  LR_inc <- vecm$model.specific$LRinclude
  lambda <- vecm$model.specific$lambda
  if(LR_inc!="none") lambda <- lambda[1:k]

  
###
  if(LR_inc=="none") {
    testCat <- switch(inc, "const"="H_lc", "both"="H_ql", "none"="H_z") 
  } else if(inc=="none"&LR_inc=="const") {
    testCat <- "H_c" 
  } else if(inc=="const"&LR_inc=="trend") {
    testCat <- "H_l"
  } else {
    stop("Sorry, rank.test does not work for model selected (due to specification of deterministic terms)\n")
  }

## trace test:
  trace <- -t*rev(cumsum(rev(log(1-lambda))))

## eigen 
  eigen <- -t*log(1-lambda)

## p-values:
  trace_pval <- gamma_doornik_all(trace, nmp=length(trace):1, test=testCat, type="trace")
  adjT <- vecm$t -floor(npar(vecm)/vecm$k)
  trace_pval_T <- gamma_doornik_all(trace, nmp=length(trace):1, test=testCat, type="trace", smallSamp=TRUE, T=adjT)
  eigen_pval <- gamma_doornik_all(eigen, nmp=length(eigen):1, test=testCat, type="eigen")

## select r
  pvals <- switch(type, trace=trace_pval, eigen=eigen_pval)
  if(missing(r_null)){
    w.pvals <- if(all(pvals<cval)) length(pvals) else which(pvals>cval)[1] 
    pval <- pvals[w.pvals]
    rank <- if(pval>cval) w.pvals-1 else w.pvals
    r_null <- "unspec"
  } else {
    pval <- pvals[r_null+1]
    rank <- r_null
    r_null <- "specified"
  }

## assemble
  res <- list()
  cal <- list(call=match.call())
  cal$cval <- cval
  cal$k <- k
  cal$type <- type
  cal$r_null <- r_null

  res_df <- data.frame(r=0:(k-1), trace=trace, trace_pval=trace_pval, trace_pval_T=trace_pval_T, eigen=eigen, eigen_pval=eigen_pval)
  res$res_df <- res_df
  res$r <- rank

  res$pval <- pval
  res$call <- cal
  class(res) <- "rank.test"
  return(res)
}

#' @rdname rank.test
#' @method print rank.test
#' @S3method print rank.test
print.rank.test <- function(x, ...) {

  if(x$call$r_null=="unspec"){
    cat("Rank selected:", x$r, "(first",x$call$type, "test with pval above", 100*x$call$cval, "%:", round(100*x$pval,1), "%)\n")
  } else {
    alter <- switch(x$call$type, "eigen" = x$r+1, "trace"= x$call$k)
    cat("Test of rank ", x$r," versus ", alter, "  (",x$call$type, " test), p-value: ", x$pval,  ".\n", sep="")
  }
  invisible(x)
}

#' @rdname rank.test
#' @param digits The number of digits to use in \code{\link{format.pval}}
#' @method summary rank.test
#' @S3method summary rank.test
summary.rank.test <- function(object, digits=max(1, getOption("digits") - 3), ...) {
  res <- object$res_df
  rownames(res) <- NULL
  res[, c(3,4,6)] <- sapply(res[, c(3,4,6)], function(x) format.pval(x, eps=1e-03,digits = digits))
  return(res)
}


#############################################
############ P val approximation
#############################################

gamma_doornik <- function(x, nmp,q,  test=c("H_z", "H_c", "H_lc", "H_l", "H_ql"), type=c("trace", "eigen"), T, smallSamp=FALSE){

  test <- match.arg(test)
  type <- match.arg(type)
  miss <- c(missing(x), missing(q))
  if(all(miss)|all(!miss)) stop("Please provide only one of args 'x' or 'q'\n")

### TRACE CASE

if(type=="trace"){
## No corr:
  paras_mean <- doornik_tab7_mean[,test]
  paras_var  <- doornik_tab7_var[,test]
  val_mean <- c(nmp^2, nmp, 1, ifelse(nmp==1,1,0), ifelse(nmp==2,1,0), sqrt(nmp))
  val_var <- c(nmp^2, nmp, 1, ifelse(nmp==1,1,0), ifelse(nmp==2,1,0))
  mean <- val_mean%*% paras_mean
  var <- val_var%*% paras_var

### Small sample case:
  if(smallSamp){
    paras_mean <- doornik_tab9_mean[,test]
    paras_var  <- doornik_tab9_var[,test]
    vals <- c(sqrt(nmp)/T,  nmp/T,  nmp^2/T^2, ifelse(nmp==1, 1/T,0), ifelse(nmp==1,1,0), ifelse(nmp==2,1,0), ifelse(nmp==3,1,0))
    mean_corr <- vals %*% paras_mean
    var_corr <- vals %*% paras_var
    mean <- exp(log(mean) +mean_corr)
    var <-  exp(log(var) +var_corr)
  }

} else {

### EIGEN CASE
  paras_mean_eig <- doornik_tab8_mean[,test]
  paras_var_eig  <- doornik_tab8_var[,test]
  val_eig <- c(nmp, 1, ifelse(nmp==1,1,0), ifelse(nmp==2,1,0), sqrt(nmp))
  mean <- val_eig%*% paras_mean_eig
  var <- val_eig%*% paras_var_eig

}

## dfs
  df1 <- mean^2/var
  df2 <- mean/var

## Compute p val or quantiles
  if(!missing(x)) res <- 1-pgamma(x, df1, df2)
  if(!missing(q)) res <- qgamma(q, df1, df2)
  names_res <- if(!missing(x)) "pval" else paste(100*q, "%", sep="")
  names(res) <- names_res

## return res
  return(res)
}


gamma_doornik_all <- function(x, nmp,q,  test=c("H_z", "H_c", "H_lc", "H_l", "H_ql"), type=c("trace", "eigen"), smallSamp=FALSE, T){

## small checks
  test <- match.arg(test)
  type <- match.arg(type)
  miss <- c(missing(x), missing(q))
  if(all(miss)|all(!miss)) stop("Provide only one of args 'x' or 'q'\n")

## many x
  if(!missing(x)){
    res <- vector("numeric", length(x))
    for(i in 1:length(x)) res[i] <- gamma_doornik(x=x[i], nmp=nmp[i], test=test, type=type, smallSamp= smallSamp, T=T)
    names(res) <-paste("pval n-p", nmp, sep="=")
  } else {
    res <- matrix(NA, ncol=length(q), nrow=length(nmp))
    for(i in 1:length(nmp)) res[i,] <- gamma_doornik(q=q, nmp=nmp[i], test=test, type=type, smallSamp= smallSamp, T=T)
    rownames(res) <- paste("n-p", nmp, sep="=")
    colnames(res) <- paste(100*q, "%", sep="")
  }
return(res)
}


### Critical values tables (taken from Doornik, in gretl plugin/johansen.c 
Hnames <- c("H_z", "H_c", "H_lc", "H_l", "H_ql")


doornik_tab7_mean <- matrix(c(2, -1, 0.07, 0.07, 0, 0, 2, 2.01, 0, 0.06, 0.05, 
0, 2, 1.05, -1.55, -0.5, -0.23, 0, 2, 4.05, 0.5, -0.23, -0.07, 
0, 2, 2.85, -5.1, -0.1, -0.06, 1.35), ncol=5, dimnames=list(1:6, Hnames))

doornik_tab7_var <- matrix(c(3, -0.33, -0.55, 0, 0,  3, 3.6, 0.75, -0.4, -0.3, 
 3, 1.8, 0, -2.8, -1.1,  3, 5.7, 3.2, -1.3, -0.5,  3, 4, 
0.8, -5.8, -2.66), ncol=5, dimnames=list(1:5, Hnames))

doornik_tab8_mean <- matrix(c(6.0019, -2.7558, 0.67185, 0.1149, -2.7764, 5.9498, 
0.43402, 0.04836, 0.018198, -2.3669, 5.8271, -1.6487, -1.6118, 
-0.25949, -1.5666, 5.8658, 2.5595, -0.34443, -0.077991, -1.7552, 
5.6364, -0.90531, -3.5166, -0.47966, -0.21447), ncol=5, dimnames=list(1:5, Hnames))

doornik_tab8_var <- matrix(c(1.8806, -15.499, 1.1136, 0.070508, 14.714, 2.2231, 
-7.9064, 0.58592, -0.034324, 12.058, 2.0785, -9.7846, -3.368, 
-0.24528, 13.074, 1.9955, -5.5428, 1.2425, 0.41949, 12.841, 2.0899, 
-5.3303, -7.1523, -0.2526, 12.393), ncol=5, dimnames=list(1:5, Hnames))

doornik_tab9_mean <- matrix(c(-0.101, 0.499, 0.896, -0.562, 0.00229, 0.00662, 0, 
0, 0.465, 0.984, -0.273, -244, 0, 0, 0.134, 0.422, 1.02, 2.17, 
-0.00182, 0, -0.00321, 0.0252, 0.448, 1.09, -0.353, 0, 0, 0, 
-0.819, 0.615, 0.896, 2.43, 0.00149, 0, 0), ncol=5, dimnames=list(1:7, Hnames))

doornik_tab9_var <- matrix(c(-0.204, 0.98, 3.11, -2.14, 0.0499, -0.0103, -0.00902, 
0.224, 0.863, 3.38, -0.807, 0, 0, -0.0091, 0.422, 0.734, 3.76, 
4.32, -0.00606, 0, -0.00718, 0, 0.836, 3.99, -1.33, -0.00298, 
-0.00139, -0.00268, -1.29, 1.01, 3.92, 4.67, 0.00484, -0.00127, 
-0.0199), ncol=5, dimnames=list(1:7, Hnames))





##########################################################################################
############ TEST: reproduce table in sec 10.1 of Doornik (1998) (with corr Doornik 1999) 
##########################################################################################

if(FALSE){
## Individual test of p-vals:
gamma_doornik(49.14, nmp=4)
gamma_doornik(19.06, nmp=3)
gamma_doornik(8.89, nmp=2)
gamma_doornik(2.35, nmp=1)

## Individual test of quantiles:
gamma_doornik(q=c(0.9, 0.95, 0.99), nmp=4)

## Vector test of p-vals/quantiles:
gamma_doornik_all(x=c(49.14, 19.06, 8.89,2.35), nmp=4:1)
gamma_doornik_all(x=c(49.14, 19.06, 8.89,2.35), nmp=4:1, test="H_lc")
gamma_doornik_all(x=c(49.14, 19.06, 8.89,2.35), nmp=4:1, test="H_l")

gamma_doornik_all(x=c(49.14, 19.06, 8.89,2.35), nmp=4:1, type="eigen")
gamma_doornik_all(x=c(49.14, 19.06, 8.89,2.35), nmp=4:1, test="H_lc", type="eigen")
gamma_doornik_all(x=c(49.14, 19.06, 8.89,2.35), nmp=4:1, test="H_l", type="eigen")

gamma_doornik_all(q=c(0.9, 0.95, 0.99), nmp=4:1)

}


##########################################################################################
############ TEST: rank.test 
##########################################################################################

if(FALSE){
library(tsDyn)
library(vars)
library(urca)

environment(rank.test) <- environment(star)

#data(Canada)
#data(denmark)

ve_can <- VECM(Canada, lag=1, estim="ML")
ve_can_l2 <- VECM(Canada, lag=2, estim="ML")

r_l1<- rank.test(ve_can)

summary(r_l1)


r_l2<- rank.test(ve_can_l2)
r_l2b<- rank.test(ve_can_l2, cval=0.15)
r_l2_ei<- rank.test(ve_can_l2, type="eigen", cval=0.25)
r_l2
r_l2b
r_l2_ei
summary(r_l2)

## with include="both"
ve_can_bo <- VECM(Canada, lag=1, estim="ML", include="both")
r_can_bo <- rank.test(ve_can_bo )
summary(r_can_bo)

## with include="both"
ve_can_none <- VECM(Canada, lag=1, estim="ML", include="none")
r_can_none <- rank.test(ve_can_none)
summary(r_can_none)

## with include="trend"
ve_can_tr <- VECM(Canada, lag=1, estim="ML", include="trend")
r_can_bo <- rank.test(ve_can_tr)

## with LRinclude=const
ve_can_LrCo <- VECM(Canada, lag=1, estim="ML", LRinclude="const")
r_can_LrCo <- rank.test(ve_can_LrCo)
summary(r_can_LrCo)

## with LRinclude=trend
ve_can_LrTr <- VECM(Canada, lag=1, estim="ML", LRinclude="trend", include="const")
ve_can_LrTr$model.specific$lambda
r_can_LrTr <- rank.test(ve_can_LrTr)
summary(r_can_LrTr )


}
