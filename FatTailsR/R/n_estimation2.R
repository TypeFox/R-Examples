

#' @include m_laplaceroll.R



#' @title Estimation and Regression Functions for Kiener Distributions
#'
#' @description
#' Several functions to estimate the parameters of asymmetric Kiener distributions 
#' and display the results in a numeric vector or in a matrix. 
#' Algorithm \code{"reg"} (the default) uses a nonlinear regression model, 
#' is slow but accurate. Algorithm \code{"estim"} just uses 5 to 11 quantiles, 
#' is very fast but less accurate.
#'
#' 
#' @param    X	       numeric. Vector, matrix, array or list of quantiles.
#' @param    algo      character. The algorithm used: \code{"r"} or \code{"reg"}  
#'                     for regression (default) and \code{"e"} or \code{"estim"}
#'                     for quantile estimation.
#' @param    ord	   integer. Option for probability selection and treatment.
#' @param    maxk	   numeric. The maximum value of tail parameter \code{k}.
#' @param    mink	   numeric. The minimum value of tail parameter \code{k}.
#' @param    maxe	   numeric. The maximum value of absolute tail parameter \code{|e|}.
#' @param    probak    numeric. Ordered vector of probabilities.
#' @param    dgts      integer. The rounding of output parameters. 
#' @param    exfitk    character. A vector of parameter names to subset the output.
#' @param    parnames  boolean. Display parameter names.
#' @param    dimnames  boolean. Display dimnames.
#' @param    ncores    integer. The number of cores for parallel processing of arrays. 
#' 
#' @details      
#' FatTailsR package currently uses two different algorithms to estimate the 
#' parameters of Kiener distributions K1, K2, K3 and K4.
#' \itemize{
#'   \item{Functions \code{fitkienerX(algo = "reg")}, \code{paramkienerX(algo = "reg")} 
#'      and \code{\link{regkienerLX}} use an unweighted  
#'      nonlinear regression from \code{logit(p)} to \code{X} over the whole dataset.  
#'      Depending the size of the dataset, calculation can be slow but is usually
#'      accurate and describes very well the last 1-10 points in the tails 
#'      (except if there is a huge outlier). }
#'   \item{Functions \code{fitkienerX(algo = "estim")}, \code{paramkienerX(algo = "estim")}, 
#'      \code{paramkienerX5} and \code{paramkienerX7} estimate the parameters with 
#'      just 5 to 11 quantiles, 5 being the minimum. For averaging purpose, 
#'      11 quantiles are proposed (see below). Computation is almost instantaneous 
#'      and reasonnably accurate. This is the recommanded method for intensive computation.}
#'   }
#' 
#' A typical input is a numeric vector or a matrix that describes the returns of a stock. 
#' A matrix must be in the format DS with DATES as rownames, STOCKS as colnames and 
#' (log-)returns as the content of the matrix. 
#' An array must be in the format DSL with DATES as rownames, STOCKS as colnames 
#' LAGS in the third dimension and (log-)returns as the content of the array. 
#' A list can be a list of numeric but neither a list of matrix, a list of data.frame 
#' or a list of arrays.
#' 
#' Conversion from a (possible) time series format to a sorted numeric vector 
#' is done automatically and without any check of the initial format. 
#' Empirical probabilities of each point in the sorted dataset is calculated 
#' with the function \code{\link{ppoints}} whose parameter \code{a} has been set to 
#' \code{a = 0} as large datasets are very common in finance. 
### The parameter \code{app} corresponds to 
### the parameter \code{a} in \code{\link[stats]{ppoints}} but has been  
### limited to range (0, 0.5). Default value is 0 as large datasets are 
### very common in finance. 
#' The lowest acceptable size of a dataset is not clear at this moment. A minimum 
#' of 11 points has been set in \code{"reg"} algorithm and a minimum of 15 points 
#' has been set in \code{"estim"} algorithm. It might change in the future. 
#' If possible, use at least 21 points. 
#' 
#' Parameter \code{algo} controls the algorithm used. Default is "reg".
#' 
#' When \code{algo = "reg"} (or \code{algo = "r"}), a nonlinear regression is performed 
#' with \code{\link[minpack.lm]{nlsLM}} from the logit of the empirical probabilities 
#' \code{logit(p)} over the quantiles X with the function \code{\link{qlkiener4}}. 
#' The maximum value of the tail parameter \code{k} is controlled by \code{maxk}.
#' An upper value \code{maxk = 10} is appropriate for datasets
#' of low and medium size, less than 20.000 or 50.000 points. For larger datasets, the
#' upper limit can be extended up to \code{maxk = 20}. When this limit is reached, 
#' the shape of the distribution is very similar to the logistic distribution 
#' (at least when \code{e = 0}) and the use of this distribution should be considered. 
#' Remember that value \code{k < 2} describes a distribution with no stable variance and 
#' \code{k < 1} describes a distribution with no stable mean.
#' 
#' When \code{algo = "estim"} (or \code{algo = "e"}),
#' 5 to 11 quantiles are used to estimate the parameters. 
#' The minimum is 5 quantiles : the median x.50, two quantiles at medium distance 
#' to the median, usually x.25 and x.75 and two quantiles located close to the extremes 
#' of the dataset, for instance x.01 and x.99 if the dataset \code{X} has more 
#' than 100 points, x.0001 and x.9999 if the dataset \code{X} has more than 
#' 10.000 points and so on if the dataset is larger. 
#' These quantiles are extracted with function \code{\link{fiveprobs}}. 
#' Small datasets must contain at least 15 different points. 
#' 
#' With the idea of averaging the results (but without any guarantee of better 
#' estimates), calculation has been extended to 11 probabilities  
#' extracted from \code{X} with the function \code{\link{elevenprobs}} where    
#' p1, p2 and p3 are the most extreme probabilities of the dataset \code{X}  
#' with values finishing either by \code{.x01} or \code{.x025} or \code{.x05}:
#' \itemize{
#'   \item{\code{p11 = c(p1, p2, p3, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p3, 1-p2, 1-p1)}}
#' }
#' 
#' Selection of subsets among these 11 probabilities is controlled with the option 
#' \code{ord} which can take 12 different values.  
#' For instance, the default \code{ord = 7} computes the  parameters at probabilities 
#' \code{c(p1, 0.25, 0.50, 0.75, 1-p1)} and \code{c(p2, 0.25, 0.50, 0.75, 1-p2)}.
#' Parameters \code{d} and \code{k} are averaged first and the results of these 
#' averages are used to compute the other parameters \code{g, a, w, e}. 
#' Small dataset should consider \code{ord = 5} and 
#' large dataset can consider \code{ord = 12}. 
#' The 12 possible values of \code{ord} are: 
#' \enumerate{
#'   \item{ \code{c(p1, 0.35, 0.50, 0.65, 1-p1)}}
#'   \item{ \code{c(p2, 0.35, 0.50, 0.65, 1-p2)}}
#'   \item{ \code{c(p1, p2, 0.35, 0.50, 0.65, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, p2, p3, 0.35, 0.50, 0.65, 1-p3, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, 0.25, 0.50, 0.75, 1-p1)}}
#'   \item{ \code{c(p2, 0.25, 0.50, 0.75, 1-p2)}}
#'   \item{ \code{c(p1, p2, 0.25, 0.50, 0.75, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, p2, p3, 0.25, 0.50, 0.75, 1-p3, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p1)}}
#'   \item{ \code{c(p2, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p2)}}
#'   \item{ \code{c(p1, p2, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p2, 1-p1)}}
#'   \item{ \code{c(p1, p2, p3, 0.25, 0.35, 0.50, 0.65, 0.75, 1-p3, 1-p2, 1-p1)}}
#' }
#'
#' \code{paramkienerX5} is a simplified version of \code{paramkienerX} with  
#' predefined values \code{algo = "estim"}, \code{ord = 5}, \code{maxk = 10} 
#' and direct access to internal subfunctions. 
#' It uses the following probabilities:
#' \itemize{
#'   \item{ \code{p5 = c(p1, 0.25, 0.50, 0.75, 1-p1)} }
#' }
#' 
#' \code{paramkienerX7} is a simplified version of \code{paramkienerX} with 
#' predefined values \code{algo = "estim"}, \code{ord = 7}, \code{maxk = 10} 
#' and direct access to internal subfunctions.
#' It uses the following probabilities:
#' \itemize{
#'   \item{ \code{p7 = c(p1, p2, 0.25, 0.50, 0.75, 1-p2, 1-p1)} }
#' }
#' 
#' The quantiles corresponding to the above probabilities are then extracted 
#' with the function \code{\link{quantile}} whose parameter \code{type} 
#' has been set to \code{type = 6} as it returns the closest values 
#' to the true quantiles (according to our experience) for all \code{k > 1.9}. 
### can change significantly the extracted quantiles. 
### Our experience is that \code{type = 6} is appropriate when \code{k > 1.9} and 
### \code{type = 5} is appropriate when \code{k < 1.9}. 
### Other types \code{type = 8} and \code{type = 9} can be considered as well. 
### The other types should be ignored. 
#' (Note: when \code{k < 1.5}, algorithm \code{algo = "reg"} returns better  
#' results). 
#' Both probabilities and quantiles are then transfered to \code{\link{estimkiener11}} 
#' for calculation.
#'  
#' \code{probak} controls the probabilities at which the model is tested with the parameter 
#' estimates. \code{fitkienerX} and \code{\link{regkienerLX}} share the same subroutines.
#' The default for \code{fitkienerX} and \code{regkienerLX} is 
#' \code{pprobs2 = c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99)} as those values 
#' are usual in finance. Other sets of values are provided at \code{\link{pprobs0}}.
#' 
#' Rounding the results is useful to display nice results, especially 
#' in a matrix or in a data.frame. \code{dgts = 13} is recommanded 
#' as \code{a}, \code{k}, \code{w} are usually significant at 1 digit.
#' \itemize{
#'   \item{ \code{dgts = NULL} does not perform any rounding. }
#'   \item{ \code{dgts = 0 to 9} rounds all parameters at the same level. }
#'   \item{ \code{dgts = 10 to 27} rounds the parameters at various levels for nice display.  
#'          See \code{\link{roundcoefk}} for the details. (Note: the
#'          rounding \code{10 to 27} currently works with \code{paramkienerX}, \code{paramkienerX5},  
#'          \code{paramkienerX7} but not yet with \code{fitkienerX}). }
#' } 
#' 
#' Extracting the most useful parameters from the (quite long) vector/matrix 
#' \code{fitk} is controlled by parameter \code{exfitk} that calls user-defined or
#' predefined parameter subsets like \code{\link{exfit0}}, ..., \code{\link{exfit7}}.
#' IMPORTANT: never subset \code{fitk} by rank number as new items may be added 
#' in the future and rank may vary.
#' 
#' Calculation of vectors, matrices and lists is not parallelized. Parallelization 
#' of code for arrays was introduced in version 1.5-0 and improved in version 1.5-1. 
#' \code{ncores} controls the number of cores allowed to the process (through 
#' \code{\link[parallel]{parApply}} which runs on Unices and Windows and requires
#' about 2 seconds to start). \code{ncores = 1} means no parallelization. 
#' \code{ncores = 0} is the recommanded option. It uses the maximum number of cores 
#' available on the computer, as detected by \code{\link[parallel]{detectCores}},  
#' minus 1 core, which gives the best performance in most cases. 
#' Although appealing, this automatic selection may be sometimes dangerous. For instance, 
#' the instruction \code{f(X, ncores_max) - f(X, ncores_max)}, a nice way to compute 
#' an array of 0, will call \code{2 ncores_max} and crash R. \code{ncores = 2,..,99} 
#' sets manually the number of cores. If the requested value is larger than the maximum 
#' number of cores, this value is automatically reduced (with a warning) to this maximum.
#' Hence, this latest option provides one core more than option \code{ncores = 0}.
#' 
#' NOTE: \code{fitkienerLX}, \code{regkienerX}, \code{estimkiener(X,5,7)} were   
#' introduced in v1.2-0 and replaced in version v1.4-1 by \code{fitkienerX} and 
#' \code{paramkiener(X,5,7)} to accomodate vector, matrix, arrays and lists. 
#' We apologize to early users who need to rewrite their codes. 
#' 
#' 
#' @return  
#' \code{paramkienerX}: a vector (or a matrix) of parameter estimates 
#' \code{c(m, g, a, k, w, d, e)}.
#' 
#' \code{fitkienerX}: a vector (or a matrix) made of several parts:
#' \itemize{
#'   \item{ \code{ret} : the return over the period calculated with \code{sum(x)}. 
#'          Thus, assume log-returns. } 
#'   \item{ \code{m, g, a, k, w, d, e} : the parameter estimates. } 
#'   \item{ \code{m1, sd, sk, ke} : the mean, standard deviation, 
#'          skewness and excess of kurtosis computed from the parameter estimates. } 
#'   \item{ \code{m1x, sdx, skx, kex} : The mean, standard deviation,  
#'          skewness and excess of kurtosis computed from the dataset. } 
#'   \item{ \code{lh} : the length of the dataset over the period. } 
#'   \item{ \code{q.} : quantile estimated with the parameter estimates. }
#'   \item{ \code{VaR.} : Value-at-Risk, positive in most cases. } 
#'   \item{ \code{c.} : corrective tail coefficient = (q - m) / (q_logistic_function - m). }
#'   \item{ \code{ltm.} : left tail mean (signed ES on the left tail, usually negative). } 
#'   \item{ \code{rtm.} : right tail mean (signed ES on the right tail, usually positive). }
#'   \item{ \code{dtmq.} : (p<=0.5 left, p>0.5 right) tail mean minus quantile. }
#'   \item{ \code{ES.} : expected shortfall, positive in most cases. }
#'   \item{ \code{h.} : corrective ES  = (ES - m) / (ES_logistic_function - m). }
#'   \item{ \code{desv.} : ES - VaR, usually positive. } 
#'   \item{ \code{l.} : quantile estimated by the tangent logistic function. }
#'   \item{ \code{dl.} : quantile - quantile_logistic_function. }
#'   \item{ \code{g.} : quantile estimated by the Laplace-Gauss function. } 
#'   \item{ \code{dg.} : quantile - quantile_Laplace_Gauss_function. }
#' }
#' 
#' IMPORTANT : if you need to subset \code{fitk}, always subset it by parameter names 
#' and never subset it by rank number as new items may be added in the future and rank may vary. 
#' Use for instance \code{\link{exfit0}}, ..., \code{\link{exfit7}}.  
#'  
#' 
#' 
#' @references
#' P. Kiener, Fat tail analysis and package FatTailsR, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' \url{http://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#
#' @seealso   \code{\link{regkienerLX}}, \code{\link{estimkiener11}}, 
#'            \code{\link{roundcoefk}}, \code{\link{exfit6}}.
#'     
#' 
#' @examples     
#' 
#' require(minpack.lm)
#' require(timeSeries)
#' 
#' ### Load the datasets and choose j in 1:16
#' DS     <- getDSdata()
#' j      <- 5
#' 
#' ### and run this block
#' probak <- c(0.01, 0.05, 0.95, 0.99)
#' X      <- DS[[j]] ; names(DS)[j]
#' elevenprobs(X)
#' fitkienerX(X, algo = "reg", dgts = 3, probak = probak)
#' fitkienerX(X, algo = "estim", ord = 5, probak = probak, dgts = 3)
#' paramkienerX(X)
#' paramkienerX5(X)
#' 
#' ### Compare the 12 values of paramkienerX(ord/row = 1:12) and paramkienerX (row 13)
#' compare <- function(ord, X) { paramkienerX(X, ord, algo = "estim", dgts = 13) }
#' rbind(t(sapply( 1:12, compare, X)), paramkienerX(X, algo = "reg", dgts = 13))
#' 
#' ### Analyze DS in one step
#' t(sapply(DS, paramkienerX, algo = "reg", dgts = 13))
#' t(sapply(DS, paramkienerX, algo = "estim", dgts = 13))
#' paramkienerX(DS, algo = "reg", dgts = 13)
#' paramkienerX(DS, algo = "estim", dgts = 13)
#' system.time(fitk_rDS <- fitkienerX(DS, algo = "r", probak = pprobs2, dgts = 3))
#' system.time(fitk_eDS <- fitkienerX(DS, algo = "e", probak = pprobs2, dgts = 3))
#' fitk_rDS
#' fitk_eDS
#' 
#' ### Subset rDS and eDS with exfit0,..,exfit7
#' fitk_rDS[,exfit4]
#' fitk_eDS[,exfit7]
#' fitkienerX(DS, algo = "e", probak = pprobs2, dgts = 3, exfitk = exfit7)
#' 
#' ### Array (new example introduced in v1.5-1)
#' ### Increase the number of cores and crash R.
#' arr <- array(rkiener1(3000), c(4,3,250))
#' paramkienerX7(arr, ncores = 2)
#' ## paramkienerX7(arr, ncores = 2) - paramkienerX(arr, ncores = 2)
#'
#' ### End
#' 
#' 
#' @export
#' @name fitkienerX
fitkienerX <- function(X, algo = c("r", "reg", "e", "estim"), ord = 7, 
                       maxk = 10, mink = 1.53, maxe = 0.5,
                       probak = pprobs2, dgts = NULL, exfitk = NULL, 
					   dimnames = FALSE, ncores = 1) {
if (maxk < 10 || maxk > 20) { 
	stop("maxk must be between 10 and 20. Can be increased with the sample size.") }
if (mink < 0.2 || mink > 2) { 
	stop("mink must be between 0.2 and 2.") }
if (ord < 1 || ord > 12) { 
	stop("ord must be between 1 and 12.") }
if (!is.element(strtrim(algo, 1)[1], c("r", "e"))) { 
	stop("algo must be r, reg, e, estim.") } 
if (!checkquantiles(probak)) { stop("probak is not ordered.") }	

cubefitkienerX <- function(X, algo, ord, maxk, mink, maxe, probak, 
                           dgts, exfitk, ncores) {
	mc <- .hnbcores(ncores)
	cl <- parallel::makeCluster(mc, methods = FALSE)
	z  <- aperm(parallel::parApply(cl, X, c(1,2), .hfitkienerX,
				algo=algo, ord=ord, type=6, 
				maxk=maxk, mink=mink, maxe=maxe, app=0, probak=probak, 
				dgts=dgts, exfitk=exfitk), c(2,1,3))
	parallel::stopCluster(cl)
return(z)
}
listfitkienerX <- function(X, algo, ord, maxk, mink, maxe, probak,  
				           dgts, exfitk, dimnames, ncores) {
	z2 <- drop(sapply(X, fitkienerX, algo, ord, maxk, mink, maxe, probak,  
				           dgts, exfitk, dimnames, ncores, simplify = "array"))
	z  <- switch(dimdimc(z2), "2" = t(z2), "3" = aperm(z2, c(3,2,1)),
	                          stop("cannot handle this format"))
return(z)
}
z <- switch(dimdimc(X),  
	"1"  = .hfitkienerX(X, algo=algo, ord=ord, type=6, 
				maxk=maxk, mink=mink, maxe=maxe, app=0, 
				probak=probak, dgts=dgts, exfitk=exfitk),
	"2"  = t(apply(X, 2, .hfitkienerX, algo=algo, ord=ord, type=6, 
				maxk=maxk, mink=mink, maxe=maxe, app=0, 
				probak=probak, dgts=dgts, exfitk=exfitk)),
	"3"  = cubefitkienerX(X, algo, ord, maxk, mink, maxe, probak,  
				          dgts, exfitk, ncores),
	"-1" = listfitkienerX(X, algo, ord, maxk, mink, maxe, probak,  
				          dgts, exfitk, dimnames, ncores),
	# "-1" = t(sapply(X, .hfitkienerX, 
				# algo=algo, ord=ord, type=6, 
				# maxk=maxk, mink=mink, maxe=maxe, app=0, 
				# probak=probak, dgts=dgts, exfitk=exfitk)),
	stop("fitkienerX cannot handle this format")
	# "3"  = aperm(apply(X, c(1,2), .hfitkienerX, 
	# "3"  = aperm(parallel::parApply(cl, X, c(1,2), .hfitkienerX,	# parallel
				# algo=algo, ord=ord, type=6, 
				# maxk=maxk, mink=mink, maxe=maxe, app=0, 
				# probak=probak, dgts=dgts, exfitk=exfitk), c(2,1,3)),
	# "-1" = t(simplify2array(parallel::mclapply(X, .hfitkienerX, 
				# algo=algo, ord=ord, type=6, 
				# maxk=maxk, mink=mink, maxe=maxe, app=0, 
				# probak=probak, dgts=dgts, exfitk=exfitk, mc.cores=mc.cores))),
	# "numeric" "matrix" "array" "list" "error"
	# "array" params x dates x stocks # aperm dates x params x stocks
	)
if (dimnames) {
	if (dimdim1(z) == 2) {
		dimnames(z) <- list("STOCKS" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]])
		}
	if (dimdim1(z) == 3) {
		dimnames(z) <- list( "DATES" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]],
							"STOCKS" = dimnames(z)[[3]]) 
		}
}
return(z)
}


.hfitkienerX <- function(X, algo = c("r", "reg", "e", "estim"), ord = 7, type = 6, 
						 maxk = 10, mink = 1.53, maxe = 0.5, app = 0, probak = pprobs2, 
						 dgts = NULL, exfitk = NULL) {

if (app < 0 || app > 0.5) { 
	stop("app (the a of ppoints) must be between 0 and 0.5. 
          Recommended values: 0, 0.375, 0.5.") }
if (maxk < 10 || maxk > 20) { 
	stop("maxk must be between 10 and 20. Can be increased with the sample size.") }
if (mink < 0.2 || mink > 2) { 
	stop("mink must be between 0.2 and 2.") }
if (ord < 1 || ord > 12) { 
	stop("ord must be between 1 and 12.") }
if (!is.element(type, c(5, 6, 8, 9))) { 
	stop("type must be 5, 6, 8 or 9.") }
if (!is.element(strtrim(algo, 1)[1], c("r", "e"))) { 
	stop("algo must be r, reg, e, estim.") }

## Data + regression or estimation => coefficients
X      <- sort(as.numeric(X[is.finite(X)])) 
coefk  <- .hparamkienerX(X, algo = algo, ord = ord, type = type, 
                         maxk = maxk, mink = mink, maxe = maxe, 
                         app = app, dgts = NULL, parnames = TRUE)
names(coefk) <- c("m","g","a","k","w","d","e") 
fitk   <- if (is.null(exfitk)) {.hfitkX(X, coefk=coefk, probak=probak, dgts=dgts)} 
                          else {.hfitkX(X, coefk=coefk, probak=probak, dgts=dgts)[exfitk]}
return(fitk)
}


.hfitkX <- function(X, coefk, probak, dgts = NULL) {

## Cumulated returns
ret            <- sum(X, na.rm = TRUE)
names(ret)     <- "ret"
## Proba pprobs2
namesk         <- getnamesk(probak)
## Moments
momk		   <- kmoments(coefk, lengthx = length(X))[c("m1", "sd", "sk", "ke")]
momx		   <- xmoments(X)[c("m1x","sdx","skx","kex","lh")]
# momx		   <- xmoments(X)[c("m1", "sd", "sk", "ke", "lh")]
# names(momx)    <-             c("m1x","sdx","skx","kex","lh")
## quantiles
quantk         <- qkiener2(p = probak, m = coefk[1], g = coefk[2], 
                                       a = coefk[3], w = coefk[5] )
names(quantk)  <- namesk$nquantk
## VaR
vark           <- varkiener2(p = probak, m = coefk[1], g = coefk[2], 
                                         a = coefk[3], w = coefk[5] )
names(vark)    <- namesk$nvark
## Coef. c.01 et consorts 
ctailk         <- ckiener2(p = probak, a = coefk[3], w = coefk[5] )
names(ctailk)  <- namesk$nctailk
## ltm
ltmk           <- ltmkiener2(p = probak, m = coefk[1], g = coefk[2], 
                                         a = coefk[3], w = coefk[5] )
names(ltmk)    <- namesk$nltmk
## rtm
rtmk           <- rtmkiener2(p = probak, m = coefk[1], g = coefk[2], 
                                         a = coefk[3], w = coefk[5] )
names(rtmk)    <- namesk$nrtmk
## dtmq (tail mean - quantile) = sign*(ES-VaR)
dtmqk          <- dtmqkiener2(p = probak, m = coefk[1], g = coefk[2], 
                                          a = coefk[3], w = coefk[5] )
names(dtmqk)   <- namesk$ndtmqk
## ES
esk            <- eskiener2(p = probak,  m = coefk[1], g = coefk[2], 
                                         a = coefk[3], w = coefk[5] )
names(esk)     <- namesk$nesk
## h.01 = ES K2 sur ES logistique
hesk           <- hkiener2(p = probak,  m = coefk[1], g = coefk[2], 
                                        a = coefk[3], w = coefk[5] )
names(hesk)    <- namesk$nhesk
## ES - VaR
desvk          <- esk - vark
names(desvk)   <- namesk$ndesvk
## logistique
logisk         <- qlogis(p = probak, location = coefk[1], scale = 2*coefk[2] )
names(logisk)  <- namesk$nlogisk
## quantile - logistique
dlogisk        <- quantk - logisk
names(dlogisk) <- namesk$ndlogisk
## Suggestion  *l pour logis => ltml, rtml, esl, dltmkl, drtmkl, deskl, h 
## => h only available in this code
## Gauss function before v1.2-50
# Xmean          <- mean(X)
# Xs             <- sd(X)
# VaR		     <- if (is.nan(Xmean) || is.nan(Xs)) {NA} else {Xmean - 2.326*Xs}
# gaussk	     <- c( VaR, -VaR +qkiener2(p=0.01, m=coefk[1], g=coefk[2], a=coefk[3], w=coefk[5])) 
# names(gaussk)  <- c("VaR", "dg.01")
## Gauss function after v1.2-50 du 27/09/2015	
gaussk	       <- qnorm(probak, mean = mean(X), sd = sd(X)) 
names(gaussk)  <- namesk$ngaussk
## quantile - Gaussienne
dgaussk	       <- quantk - gaussk 
names(dgaussk) <- namesk$ndgaussk 
## Final object fitk + arrondi
fitk           <- c(ret, coefk, momk, momx, quantk, vark, ctailk, ltmk, rtmk, dtmqk, 
                    esk, hesk, desvk, logisk, dlogisk, gaussk, dgaussk)
if (!is.null(dgts)) { fitk <- round(fitk, digits = dgts) }
return(fitk)
}
 
#' @export
#' @rdname fitkienerX
paramkienerX <- function(X, algo = c("r", "reg", "e", "estim"), 
                         ord = 7, maxk = 10, mink = 1.53, maxe = 0.5, dgts = 3, 
						 parnames = TRUE, dimnames = FALSE, ncores = 1) {
if (maxk < 10 || maxk > 20) { 
	stop("maxk must be between 10 and 20. Can be increased with the sample size.") }
if (mink < 0.2 || mink > 2) { 
	stop("mink must be between 0.2 and 2.") }
if (ord < 1 || ord > 12) { 
	stop("ord must be between 1 and 12.") }
if (!is.element(strtrim(algo, 1)[1], c("r", "e"))) { 
	stop("algo must be r, reg, e, estim.") }

cubekienerX <- function(X, algo, ord, maxk, mink, maxe, dgts, 
                        parnames, ncores) {
	mc <- .hnbcores(ncores)
	cl <- parallel::makeCluster(mc, methods = FALSE)
	z  <- aperm(parallel::parApply(cl, X, c(1,2), .hparamkienerX,
					algo=algo, ord=ord, type=6, 
					maxk=maxk, mink=mink, maxe=maxe, app=0, 
					dgts=dgts, parnames=parnames), c(2,1,3))
	parallel::stopCluster(cl)
return(z)
}
listkienerX <- function(X, algo, ord, maxk, mink, maxe, dgts, 
						parnames, dimnames, ncores) {
	z2 <- drop(sapply(X, paramkienerX, algo, ord, maxk, mink, maxe, dgts, 
				  parnames, dimnames, ncores, simplify = "array"))
	z  <- switch(dimdimc(z2), "2" = t(z2), "3" = aperm(z2, c(3,2,1)),
	                          stop("cannot handle this format"))
return(z)
}
z <- switch(dimdimc(X),  
	"1" = .hparamkienerX(X, algo=algo, ord=ord,
					maxk=maxk, mink=mink, maxe=maxe, dgts=dgts, 
					parnames=parnames),
	"2"  = t(apply(X, 2, .hparamkienerX, algo=algo, ord=ord, 
					maxk=maxk, mink=mink, maxe=maxe, dgts=dgts, 
					parnames=parnames)),
	"3"  = cubekienerX(X, algo, ord, maxk, mink, maxe, dgts, 
					parnames, ncores),
	"-1" = listkienerX(X, algo, ord, maxk, mink, maxe, dgts, 
					parnames, dimnames, ncores),
	stop("paramkienerX cannot handle this format")
	# "3"   = aperm(apply(X, c(1,2), .hparamkienerX, 
					# algo=algo, ord=ord, type=6, 
					# maxk=maxk, mink=mink, maxe=maxe, app=0, dgts=dgts, 
					# parnames=parnames), c(2,1,3)),
	# "3"  = aperm(parallel::parApply(cl, X, c(1,2), .hparamkienerX,	# parallel
					# algo=algo, ord=ord, 
					# maxk=maxk, mink=mink, maxe=maxe, dgts=dgts, 
					# parnames=parnames), c(2,1,3)),
	# "-1" = t(simplify2array(parallel::mclapply(X, .hparamkienerX, 
	                # algo=algo, ord=ord, type=6, 
					# maxk=maxk, mink=mink, maxe=maxe, app=0, dgts=dgts, 
					# parnames=parnames, mc.cores=mc.cores))),
	# "-1" = t(sapply(X, .hparamkienerX, 
					# algo=algo, ord=ord, 
					# maxk=maxk, mink=mink, maxe=maxe, dgts=dgts, 
					# parnames=parnames)),
	# "numeric" "matrix" "array" "list" "error"
	# "array" params x dates x stocks # aperm dates x params x stocks
	)
if (dimnames) {
	if (dimdim1(z) == 2) {
		dimnames(z) <- list("STOCKS" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]])
		}
	if (dimdim1(z) == 3) {
		dimnames(z) <- list( "DATES" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]],
							"STOCKS" = dimnames(z)[[3]]) 
		}
}
return(z)
}


.hparamkienerX <- function(X, algo = c("r", "reg", "e", "estim"), ord = 7, type = 6, 
           maxk = 10, mink = 1.53, maxe = 0.5, app = 0, dgts = NULL, parnames = TRUE) {
		   
# library("minpack.lm") # ajoute pour parallelisation
X    <- sort(as.numeric(X[is.finite(X)]))

if (strtrim(algo, 1)[1] == "e") { 
	if (length(X) > 14) { 
			p11   <- elevenprobs(X)
			x11   <- quantile(X, probs = p11, na.rm = TRUE, names = FALSE, type = type)
			coefk <- estimkiener11(x11, p11, ord, maxk) 
		} else { 
			coefk <- rep(NA, 7) 
	}
} else {  # behaviour: (strtrim(algo, 1)[1] != "e" rather than == "r")
	if (length(unique(X)) > 10) {
		L      <- logit(ppoints(length(X), a = app)) # a = app 
		dfrXL  <- data.frame(X=X, L=L)
		## Initialisation
		parini <- .hparamkienerX5(X, parnames = FALSE)
		if (anyNA(parini)) { 
				mini   <- median(X)
				gini   <- 0.25*sd(X)
				qqq    <- quantile(X, c(0.10, 0.50, 0.90), type = 6)
				dini   <- if (anyNA(qqq)) {0} else {log(abs(qqq[3]-qqq[2])/abs(qqq[2]-qqq[1]))/4.394}
				kini   <- 4
				eini   <- min(max(-maxe, dini*kini), maxe)
				aini   <- kini/(1-eini)
				wini   <- kini/(1+eini)
			} else {
				mini   <- parini[1]
				gini   <- parini[2]
				aini   <- parini[3]
				kini   <- min(max(mink, parini[4]), maxk)
				wini   <- parini[5]
				dini   <- parini[6]
				eini   <- min(max(-maxe, parini[7]), maxe)
			}
		## Regression K4
		regk0  <- minpack.lm::nlsLM( X ~ FatTailsR::qlkiener4(L, mini, g, k, e), 
		# regk0  <- nlsLM( X ~ qlkiener4(L, mini, g, k, e), 
						 data = dfrXL, 
						 start = list(g = gini, k = kini, e = eini), 
						 lower = c(gmin = 0,   kmin = mink, emin =-maxe), 
						 upper = c(gmax = Inf, kmax = maxk, emax = maxe) 
						)
		coefk  <-    c(m = mini,
					   g = coef(regk0)[1],
					   a = ke2a(coef(regk0)[2], coef(regk0)[3]),
					   k = coef(regk0)[2],
					   w = ke2w(coef(regk0)[2], coef(regk0)[3]),
					   d = ke2d(coef(regk0)[2], coef(regk0)[3]),
					   e = coef(regk0)[3]
						) 
		## Regression K2 (maybe one day) or K1
		# regk0  <- minpack.lm::nlsLM( X ~ qlkiener2(L, mini, g, a, w), 
						 # data = dfrXL, 
						 # start = list(g = gini,   a = aini,    w = wini), 
						 # lower = c(gmin = 0,   amin = mink, wmin = mink), 
						 # upper = c(gmax = Inf, amax = maxk, wmax = maxk) 
						# )
		# coefk  <-    c(m = mini,
					   # g = coef(regk0)[1],
					   # a = coef(regk0)[2],
					   # k = aw2k(coef(regk0)[2], coef(regk0)[3]),
					   # w = coef(regk0)[3],
					   # d = aw2d(coef(regk0)[2], coef(regk0)[3]),
					   # e = aw2e(coef(regk0)[2], coef(regk0)[3])
						# ) 
	} else {
		coefk  <- rep(NA, 7)		
	}
}
z  <- roundcoefk(coefk, dgts, parnames)
return(z)
}


#' @export
#' @rdname fitkienerX
paramkienerX7 <- function(X, dgts = 3, parnames = TRUE, dimnames = FALSE, ncores = 1) { 
cubekienerX7 <- function(X, dgts, parnames, ncores) {
	mc <- .hnbcores(ncores)
	cl <- parallel::makeCluster(mc, methods = FALSE)
	z  <- aperm(parallel::parApply(cl, X, c(1,2), 
	            .hparamkienerX7, dgts, parnames), c(2,1,3))
	parallel::stopCluster(cl)
return(z)
}
listkienerX7 <- function(X, dgts, parnames, dimnames, ncores) {
	z2 <- drop(sapply(X, paramkienerX7, dgts, parnames, 
	                  dimnames, ncores, simplify = "array"))
	z  <- switch(dimdimc(z2), "2" = t(z2), "3" = aperm(z2, c(3,2,1)),
	                          stop("cannot handle this format"))
return(z)
}
z <- switch(dimdimc(X),  
	 "1" = .hparamkienerX7(X, dgts, parnames),
	 "2" = t(apply(X, 2, .hparamkienerX7, dgts, parnames)),
	 "3" = cubekienerX7(X, dgts, parnames, ncores),
	"-1" = listkienerX7(X, dgts, parnames, dimnames, ncores),
	stop("paramkienerX7 cannot handle this format")
	 # "3" = aperm(apply(X, c(1,2), .hparamkienerX7, dgts, parnames), c(2,1,3)),
	# "-1" = t(simplify2array(parallel::mclapply(X, .hparamkienerX7, 
				# dgts, parnames, mc.cores=mc.cores))),
	# "-1" = t(sapply(X, .hparamkienerX7, dgts, parnames)),
	# "numeric" "matrix" "array" "list" "error"
	# "array" params x dates x stocks # aperm dates x params x stocks
	)
if (dimnames) {
	if (dimdim1(z) == 2) {
		dimnames(z) <- list("STOCKS" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]])
		}
	if (dimdim1(z) == 3) {
		dimnames(z) <- list( "DATES" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]],
							"STOCKS" = dimnames(z)[[3]]) 
		}
}
return(z)
}


.hparamkienerX7 <- function(X, dgts = NULL, parnames = TRUE) { 
	X   <- sort(as.numeric(X[is.finite(X)]))
	if (length(X) > 14) { 
		p7    <- sevenprobs(X)
		x7    <- quantile(X, probs = p7, na.rm = TRUE, names = FALSE, type = 6) 
		coefk <- estimkiener7(x7, p7) 
	} else { 
		coefk <- rep(NA, 7) 
	}
z  <- roundcoefk(coefk, dgts, parnames)
return(z)
} 


#' @export
#' @rdname fitkienerX
paramkienerX5 <- function(X, dgts = 3, parnames = TRUE, dimnames = FALSE, ncores = 1) { 
cubekienerX5 <- function(X, dgts, parnames, ncores) {
	mc <- .hnbcores(ncores)
	cl <- parallel::makeCluster(mc, methods = FALSE)
	z  <- aperm(parallel::parApply(cl, X, c(1,2), 
	            .hparamkienerX5, dgts, parnames), c(2,1,3))
	parallel::stopCluster(cl)
return(z)
}
listkienerX5 <- function(X, dgts, parnames, dimnames, ncores) {
	z2 <- drop(sapply(X, paramkienerX5, dgts, parnames, 
	                  dimnames, ncores, simplify = "array"))
	z  <- switch(dimdimc(z2), "2" = t(z2), "3" = aperm(z2, c(3,2,1)),
	                          stop("cannot handle this format"))
return(z)
}
z <- switch(dimdimc(X) ,  
	 "1" = .hparamkienerX5(X, dgts, parnames),
	 "2" = t(apply(X, 2, .hparamkienerX5, dgts, parnames)),
	 "3" = cubekienerX5(X, dgts, parnames, ncores),
	"-1" = listkienerX5(X, dgts, parnames, dimnames, ncores),
	stop("paramkienerX5 cannot handle this format")
	 # "3" = aperm(apply(X, c(1,2), .hparamkienerX5, dgts, parnames), c(2,1,3)),
	# mc.cores <- if(tolower(.Platform$OS.type) == "windows") {1} else {mc}
	# "-1" = t(simplify2array(parallel::mclapply(X, .hparamkienerX5, 
				# dgts, parnames, mc.cores=mc.cores))),
	# "-1" = t(sapply(X, .hparamkienerX5, dgts, parnames)),
	# "numeric" "matrix" "array" "list" "error"
	# "array" params x dates x stocks # aperm dates x params x stocks
	)
if (dimnames) {
	if (dimdim1(z) == 2) {
		dimnames(z) <- list("STOCKS" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]])
		}
	if (dimdim1(z) == 3) {
		dimnames(z) <- list( "DATES" = dimnames(z)[[1]],
							"PARAMS" = dimnames(z)[[2]],
							"STOCKS" = dimnames(z)[[3]]) 
		}
}
return(z)
}


.hparamkienerX5 <- function(X, dgts = NULL, parnames = TRUE) { 
	X   <- sort(as.numeric(X[is.finite(X)]))
	if (length(X) > 14) { 
		p5    <- fiveprobs(X)
		x5    <- quantile(X, probs = p5, na.rm = TRUE, names = FALSE, type = 6) 
		coefk <- estimkiener5(x5, p5) 
	} else {
		coefk <- rep(NA, 7) 
	}
z  <- roundcoefk(coefk, dgts, parnames)
return(z)
}


.hnbcores  <- function(ncores = 1) {
	if (is.null(ncores)) {
		warning("NULL cores requested. Reverts to 1 core.")
		ncores <- 1
	}
	ncores <- ncores[1]
	if (!is.element(ncores, 0:99)) {
		warning("ncores poorly defined and not in 0:99. Reverts to 1 core.")
		ncores <- 1
	}
	if (ncores == 1) {z <- 1} else {
		ncmax <- parallel::detectCores()
		if (is.na(ncmax))   {warning("NA cores detected. Reverts to 1 core.")
							 ncmax <- 1}
		if (ncores > ncmax) {warning(
		 paste0(sprintf("%d", ncores), 
				" cores requested. Reverts to the maximum of this processor: ", 
				sprintf("%d", ncmax), 
				" cores."))
		 }
		if (ncores == 0) {z <- max(1, ncmax - 1)} else {z <- min(ncores, ncmax)}
	}
return(z)
}

