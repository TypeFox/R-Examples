#-----------------------------------------------------------------------------#
#                                                                             #
#                     QUALITY CONTROL STATISTICS IN R                         #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores S?nchez                                       #
#              Student Master of Statistical Techniques                       #
#              University of The Coru?a, SPAIN                                #
#              mflores@outlook.com                                            #
#                                                                             #
#-----------------------------------------------------------------------------#
#-------------------------------------------------------------------------
# mewma chart
#-------------------------------------------------------------------------
##' Function to plot mewma chart
##'
##' This function is used to compute statistics required by the mewma chart.
##'
##' @param x   an R object (used to select the method). See details.
##' @param ... arguments passed to or from methods.
##' @export
## @references Montgomery, D.C. (2000)
##' @examples
##' 
##' ##
##' ##  Continuous data 
##' ##
##' library(qcr)
##' data(dowel1)
##' str(dowel1)
##' data.mqcd <- mqcd(dowel1)
##' res.mqcs <- mqcs.mewma(data.mqcd)
##' summary(res.mqcs)
##' plot(res.mqcs, title =" MEWMA Control Chart for dowel1")

mqcs.mewma <- function(x, ...) {
  UseMethod("mqcs.mewma")
}

##' @rdname mqcs.mewma
##' @method mqcs.mewma default
##' @inheritParams mqcd
##' @param Xmv is the mean vector. It is only specified for Phase II or when the parameters of the distribution are known.
##' @param S is the sample covariance matrix. It is only used for Phase II or when the parameters of the distribution are known.
##' @param lambda is the smoothing constant. Only values of 0.1, 0.2,...,0.9 are allowed.
##' @param method Is the method employed to compute the covatiance matrix
##' in individual observation case. Two methods are used "sw" 
##' for compute according to (Sullivan,Woodall 1996a) and "hm" 
##' by (Holmes,Mergen 1993)
##' @param plot a logical value indicating should be plotted. 
##' @author Edgar Santos-Fernandez
##' @export
##' 
mqcs.mewma.default <- function(x, data.name = NULL, Xmv = NULL, S = NULL,
                            method = "sw", plot = FALSE, ...)
#.........................................................................
  {
  
  obj<-mqcd(data= x, data.name = data.name)

  result<-mqcs.mewma.mqcd(x = obj, 
                          Xmv = Xmv, S = S, lambda = 0.1,
                          method = method, plot = plot, ...)

  return(result)
} # mqcs.mewma.default
#.........................................................................

##' @rdname  mqcs.mewma
##' @method mqcs.mewma mqcd
##' @inheritParams mqcs.mewma.default
##' @export
##' 

mqcs.mewma.mqcd <- function(x, Xmv = NULL, S = NULL, 
                            lambda = 0.1,
                            method = "sw", plot = FALSE, ...) 
#.........................................................................  
{
  
  if(is.null(x) || !inherits(x, "mqcd"))
    stop("data must be an objects of class (or extending) 'mqcd'")
    
  mqcs<-mqcs(x, method)
  if(is.null(Xmv)) Xmv <- mqcs$mean 
  if(is.null(S)) S <- mqcs$S
  x.jk <- mqcs$mean.jk
  
  p <- ncol(x) # quality characteristics
  m <- nrow(x) # number of samples or observations
  n <- dim(x)[3] # observations or sample size 
  
  
  statistics <- matrix(0,m,1)
  
  h4 <- matrix(c(8.6336,9.6476,10.083,10.3114,10.4405,10.5152,10.5581,10.5816,10.5932,10.814,
                 11.8961,12.3505,12.5845,12.7143,12.788,12.8297,12.8524,12.8635,12.7231,13.8641,14.3359,
                 14.576,14.7077,14.7818,14.8234,14.846,14.857,14.5363,15.7293,16.217,16.4629,16.5965,16.6711,
                 16.7127,16.7352,16.7463,16.2634,17.5038,18.0063,18.2578,18.3935,18.4687,18.5105,18.5331,
                 18.5442,17.9269,19.2113,19.7276,19.9845,20.1223,20.1982,20.2403,20.2631,20.2743,19.541,
                 20.8665,21.396,21.6581,21.798,21.8747,21.9171,21.9401,21.9515,21.1152,22.4796,23.0217,23.2887,
                 23.4307,23.5082,23.551,23.5742,23.5858,22.6565,24.0579,24.6119,24.8838,25.0278,25.1062,25.1493,
                 25.1728,25.1846),nrow = 9)
  
  rownames(h4) <- c(seq(0.1,0.9, by  = 0.1))
  colnames(h4) <- c(1:9)
  
  z<-matrix(0, m, p)
  m1 <- rownames(h4)
  m2 <- colnames(h4)
  l <- lambda * 10
  ucl <- h4[m1[l], m2[p - 1]]

  
  for (i in 1 : m){
    if(i==1){
      z[i,] <- lambda * (x.jk[i,] - Xmv)}
    else{ z[i,] <- lambda * (x.jk[i,] - Xmv) + (1 - lambda) * z[i - 1,]}
    weig <- S * (lambda * (1 - ((1 - lambda) ^ (2 * i))) / (2 - lambda))
    statistics[i,1] <- t(z[i,]) %*% solve(weig) %*% z[i,]
  }
  
  limits <- c(lcl = 0, ucl = ucl)
  
  violations <- which(statistics > ucl)
    

  data.name <- attr(x, "data.name")
  result <- list(mqcd  =  x, type  =  "mewma", statistics  =  statistics,
                 mean  =  Xmv, S  =  S, lambda = lambda,
                 limits  =  limits, data.name  =  data.name,
                 violations  =  violations)
  
  oldClass(result) <- c("mqcs.mewma", "mqcs")
  
  if(plot) plot(result, ...)
  
  return(result)
#.........................................................................
} # mqcs.mewma.mqcd
#.........................................................................