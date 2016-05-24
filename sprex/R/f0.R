#' @name f0
#' @title Number of Unobserved Species
#' @description Calculate the number of unobserved species (f0).
#' 
#' @param f a vector of species frequencies where \code{f[i]} is the number of 
#'   species represented by only \code{i} samples.
#' @param N population size.
#'  
#' @return All functions return a vector containing the estimated number of 
#'   species (\code{s.est}), unobserved species (\code{f0}), observed species 
#'   (\code{s.obs}), and the total number of samples (\code{n}). \code{Swor1} 
#'   also returns the standard deviation of \code{s.est} as \code{sd.s.est}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @references 
#'   \code{Chao1,ACE}: Colwell, R.K., A. Chao, N.J. Gotelli, S.-Y. Lin, 
#'   C.X. Mao, R.L. Chazdon, and J.T. Longino. 2012. Models and estimators 
#'   linking individual-based and sample-based rarefaction, extrapolation and 
#'   comparison of assemblages. Journal of Plant Ecology 5(1):3-21. \cr\cr
#'   \code{jack1,jack2}: Burnham, KP and WS Overton. 1978. Estimation of the 
#'   size of a closed population when capture probabilities vary among animals. 
#'   Biometrika 65(3):625-633. \cr\cr
#'   \code{Swor1}: Chao, A. and C.-W. Lin. 2012. Nonparametric lower bounds for 
#'   species richness and shared species richness under sampling without 
#'   replacement. Biometrics 68:912-921. \cr\cr
#'   \code{iChao1}: Chiu, C-H, Wang, Y-T, Walther, BA, and A Chao. 2014. 
#'   An impro.ved nonparametric lower bound of species richness via a 
#'   modified Good-Turing frequency formula. Biometrics 70(3):671-682. \cr\cr
#'   \code{clench}: Clench, H. 1979. How to make regional lists of butterflies: 
#'   Some thoughts. Journal of the Lepidopterists' Society 33(4):216-231
#'   
#' @examples
#' data(osa.second.growth)
#' f <- expand.freqs(osa.second.growth)
#' 
#' ace.est <- ACE(f)
#' chao1.est <- Chao1(f)
#' jack1.est <- jack1(f)
#' jack2.est <- jack2(f)
#' swor1.est <- Swor1(f, 20000)
#' ichao1.est <- iChao1(f)
#' clench.est <- Clench(f, num.reps = 50)
#' 
#' f0.est <- cbind(
#'   ACE = ace.est["f0"],
#'   Chao1 = chao1.est["f0"],
#'   jack1 = jack1.est["f0"],
#'   jack2 = jack2.est["f0"],
#'   Swor1 = swor1.est["f0"],
#'   iChao1 = ichao1.est["f0"],
#'   clench = clench.est["f0"]
#' )
#' f0.est
#' 
#' @aliases f0
#' 
NULL