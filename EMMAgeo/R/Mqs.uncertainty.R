#' Function to estimate end-member scores uncertainty
#' 
#' The function uses either existing assemblages of end-member loadings or
#' specified measures of centrality and dispersion as input for Monte Carlo
#' runs to estimate the influence of different end-member loadings on
#' end-member scores. Likewise, the influence of the weight limit quantiles
#' (l) can be estimated.
#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param q Numeric scalar with the number of end-members to include. Only
#' necessary in combination with \code{Vqn} as matrix of user-defined
#' end-member loadings.
#' @param l Numeric vector with the weight tranformation limits (i.e.
#' quantiles after Klovan & Imbrie, 1971). If the parameter is of length 1,
#' \code{l} is assumed to be a constant, if of length 2, \code{l} defines
#' either mean and standard deviation or minimum and maximum, depending on the
#' value of \code{type.l}.
#' @param c Numeric scalar specifying the constant sum scaling parameter, e.g.
#' 1, 100, 1000; default is 100.
#' @param rotation Character scalar, rotation type, default is "Varimax" (cf.
#' Dietze et al., 2012). One out of the rotations provided in GPArotation is
#' possible (cf. \code{\link{rotations}}).
#' @param Vqn Numeric matrix with existing unscaled end-member loadings. These
#' may represent user-defined loadings (or mean loadings if \code{Vqn.sd} is
#' specified). See example section for further information.
#' @param Vqn.sd Numeric matrix with standard deviations of the mean unscaled
#' end-member loadings in \code{Vqn}.
#' @param runs Logical scalar with the number of Monte Carlo runs to be
#' performed, default is 100.
#' @param type.l Character scalar with type of random l value generation.
#' Either \code{"rnorm"} or \code{"runif"}, default is \code{"runif"}.
#' @param autocorrelation Numeric scalar optionally specifying the degree of
#' autocorrelation among classes. Autocorrelation is realised as running mean
#' of the specified length. Only odd values are allowed.
#' @param \dots Further arguments passed to the function.
#' @return A list with numeric vector and matrix objects. \item{l}{Randomised
#' weight limit values.} \item{Vqn}{Randomised unscaled end-member loadings.}
#' \item{Mqs}{Modelled end-member scores.} \item{mean}{Modelled end-member
#' score means.} \item{sd}{Modelled end-member score standard deviations.}
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{test.robustness}}, \code{\link{test.parameters}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.\cr
#' Klovan JE, Imbrie J. 1971. An Algorithm and FORTRAN-IV Program for
#' Large-Scale Q-Mode Factor Analysis and Calculation of Factor Scores.
#' Mathematical Geology 3: 61-77.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X, envir = environment())
#' 
#' ## set model run parameters
#' q = 3 # set number of end-members, try 4 to see the difference!
#' Vqn <- EMMA(X, q)$Vqn # assign unscaled end-member loadings
#' Vqn.sd <- Vqn * 0.2 # assign a relative standard deviation of 20 %
#' l.1 <- 0.2 # set l to 0.2
#' l.2 <- c(0.2, 0.08) # set l to mean = 0.2 and sd = 0.08
#' runs <-  12 # senseless value to increase computation speed
#' 
#' ## EXAMPLE 1
#' ## Calculate Mqs uncertainty
#' M <- Mqs.uncertainty(X = X, 
#'                      q = q, 
#'                      l = l.1,
#'                      runs = runs,
#'                      Vqn = Vqn,
#'                      Vqn.sd = Vqn.sd,
#'                      type.l = "rnorm",
#'                      autocorrelation = 3)
#' 
#' ## Plot line-point graph with means and standard deviations
#' plot(NA,
#'      xlim = c(1, nrow(X)),
#'      ylim = c(0.5, q + 1),
#'      main = "End-member scores with uncertainty")
#' for(i in 1:q) {
#'   lines(1:nrow(X), M$mean[,i] - M$sd[,i] + i, col = i, lty = 2)
#'   lines(1:nrow(X), M$mean[,i] + i, col = i, lwd = 2)
#'   points(1:nrow(X), M$mean[,i] + i, col = i)
#'   lines(1:nrow(X), M$mean[,i] + M$sd[,i] + i, col = i, lty = 2)
#' }
#' 
#' ## EXAMPLE 2
#' ## Calculate Mqs uncertainty
#' M <- Mqs.uncertainty(X = X, 
#'                      q = q, 
#'                      l = l.2,
#'                      runs = runs,
#'                      Vqn = Vqn,
#'                      type.l = "rnorm")
#' 
#' ## Plot point graph with error bars
#' plot(NA,
#'      xlim = c(1, nrow(X)),
#'      ylim = c(0.5, q + 1),
#'      main = "End-member scores with uncertainty")
#' for(i in 1:q) {
#'   points(1:nrow(X), M$mean[,i] + i, pch = 3, col = i)
#'   arrows(1:nrow(X), M$mean[,i] - M$sd[,i] + i, 
#'          1:nrow(X), M$mean[,i] + M$sd[,i] + i, 
#'          code = 3, angle = 90, length = 0.05, col = i)
#' }
#' 
#' @export Mqs.uncertainty
Mqs.uncertainty <- function(
  X,
  q,
  l,
  c,
  rotation = "Varimax",
  Vqn,
  Vqn.sd,
  runs,
  type.l,
  autocorrelation,
  ...
){

  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  if("type.lw" %in% names(list(...))) {
    stop('Parameter "type.lw" is depreciated. Use "type.l" instead.')
  }
  
  ## test data consistency
  if(missing(Vqn) == FALSE &
    (missing(q) == TRUE |
     missing(l) == TRUE)) {stop(paste("Vqn.data present but",
       "no value for q or l is provided. Cannot execute opration!"))}
  
  ## check/set default values
  if(missing(c) == TRUE) {c <- 100}
  if(missing(runs) == TRUE) {runs <- 1000}
  if(missing(type.l) == TRUE) {type.l = "runif"}
  
  ## define output variables
  Mqs.out <- array(dim = c(nrow(X), q, runs))
  Mqs.mean <- matrix(nrow = nrow(X), ncol = q)
  Mqs.sd <- matrix(nrow = nrow(X), ncol = q)
  
  ## generate random data
  ## l data set
  if(length(l) == 1) {
    l.out <- rep(l, runs)
  } else if(length(l) == 2) {
    if(type.l == "runif") {
      l.out <- runif(n = runs, min = l[1], max = l[2])
    } else if(type.l == "rnorm") {
      l.out <- rnorm(n = runs, mean = l[1], sd = l[2])
    }
  }
  l.out <- ifelse(l.out < 0, 0, l.out) # set negative values to zero
  
  ## Vqn data set
  if(missing(Vqn.sd) == TRUE) {
    Vqn.out <- rep(Vqn, runs)
    dim(Vqn.out) <- c(q, ncol(X), runs)
  } else {
    Vqn.out <- array(dim = c(q, ncol(X), runs))
    for(i in 1:runs) {
      for(j in 1:q) {
        Vqn.out[j,,i] <- rnorm(n = ncol(X), 
                              mean = Vqn[j,],
                              sd = Vqn.sd[j,])
        if(missing(autocorrelation) != TRUE) {
          if(autocorrelation %% 2 == 0) {
            stop("Value for autocorrelation is no odd integer.")
          }
          b.h <- autocorrelation / 2 - 0.5
          running.mean <- Vqn.out[j,,i]
          for(k in (1 + b.h):(ncol(X) - b.h)) {
            running.mean[k] <- mean(Vqn.out[j,(k - b.h):(k + b.h),i], 
                                    na.rm = TRUE)
          }
          Vqn.out[j,,i] <- running.mean
        }
      }
    }
  }
  
  ## perform EMMA with specified input data
  for(i in 1:runs) {
    Mqs.out[,,i] <- EMMA(X = X,
                         q = q,
                         l = l.out[i],
                         c = c,
                         Vqn = Vqn.out[,,i],
                         rotation = rotation
                         )$Mqs
  }
  
  ## calculate statistical summary
  for(i in 1:q) {
    Mqs.mean[,i] <- apply(Mqs.out[,i,], 1, mean, na.rm = TRUE)
    Mqs.sd[,i] <- apply(Mqs.out[,i,], 1, sd, na.rm = TRUE)
  }
  
  ## rescale Mqs.mean to 100 %
  Mqs.mean <- Mqs.mean / apply(Mqs.mean, 1, sum, na.rm = TRUE)

  return(list(l = l.out,
              Vqn = Vqn.out,
              Mqs = Mqs.out,
              mean = Mqs.mean,
              sd = Mqs.sd))
}