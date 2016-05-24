

#' @include i_kiener4.R



#' @title Moments Associated To Kiener Distribution Parameters
#'
#' @description
#' Non-central moments, central moments, mean, standard deviation, variance, 
#' skewness, kurtosis, excess of kurtosis and cumulants associated to  
#' the parameters of Kiener distributions K1, K2, K3 and K4. 
#' All-in-one vectors \code{kmoments} (estimated from the parameters) 
#' and \code{xmoments} (estimated from the vector of quantiles) are provided.
#' 
#' @param    x	    	numeric. Vector of quantiles.
#' @param    n	    	integer. The moment order. 
#' @param    coefk    	vector. Parameters of the distribution of length 3 ("K1"),
#'                      length 4 (model = K2, K3, K4) and length 7 ("K7"). 
#' @param    model   	character. Model type, either "K2", "K3" or "K4" if \code{coefk} is 
#'                      of length 4. Type "K1" and "K7" may be provided but are ignored. 
#' @param    lengthx	integer. The length of the vector \code{x} used to calculate the parameters.
#'                      See the details for matrix and lists.
#' @param    dgts       integer. The rounding applied to the output.
#' @param    dimnames   boolean. Display dimnames.
#' 
#' @details 
#' The non-central moments \code{m1,m2,m3,m4,..,mn}, 
#' the central moments \code{u1,u2,u3,u4,..,un} (where u stands for mu in Greek)
#' and the cumulants \code{k1,k2,k3,k4,..,kn} (where k stands for kappa in Greek; 
#' not to be confounded with tail parameter "k" and models "K1", "K2", "K3", "K4") 
#' of order \eqn{n} exist only if \eqn{min(a, k, w) > n}. 
#' The mean \code{m1} exists only if \eqn{min(a, k, w) > 1}. 
#' The standard deviation \code{sd} and the variance \code{u2} exist only
#' if \eqn{min(a, k, w) > 2}. 
#' The skewness \code{sk} exists only if \eqn{min(a, k, w) > 3}. 
#' The kurtosis \code{ku} and the excess of kurtosis \code{ke} exist only 
#' if \eqn{min(a, k, w) > 4}. 
#' 
#' \code{coefk} may take five different forms :
#' \itemize{
#'   \item{\code{c(m, g, k) } of length 3 for distribution "K1".}
#'   \item{\code{c(m, g, a, w) } of length 4 for distribution "K2".}
#'   \item{\code{c(m, g, k, d) } of length 4 for distribution "K3".}
#'   \item{\code{c(m, g, k, e) } of length 4 for distribution "K4".}
#'   \item{\code{c(m, g, a, k, w, d, e) } of length 7 (sometimes referred as "K7") 
#'         provided by estimation/regression functions \code{paramkienerX}, 
#'         \code{fitkienerX}, \code{regkienerLX} (via \code{"reg$coefk"})
#'         and conversion function \code{pk2pk}.}        
#' }
#' Forms of length 3 and 7 are automatically recognized and do not require 
#' \code{model = "K1"} or \code{"K7"} which are ignored. 
#' Forms of length 4 require \code{model = "K2"}, \code{"K3"} or \code{"K4"}. 
#' Visit \code{\link{pk2pk}} for details on the parameter conversion function 
#' used within \code{kmoments}.
#' 
#' \code{xmoments} and \code{kmoments} provide all-in-one vectors. 
#' 
#' \code{xmoments} is the traditional mean of squares, cubic and power 4 functions 
#' of non-central and central values of x, from which NA values have been removed. 
#' Therefore, length of x ignores NA values and may be different from the true length.
#' 
#' \code{kmoments} calls every specialized functions from order 1 to order 4 and 
#' uses the estimated parameters as inputs, not the initial dataset \code{x}.
#' As it does not know \emph{a priori} the length of \code{x}, this latest can 
#' be provided separately via \code{lengthx = length(x)}, \code{lengthx = nrow(x)} 
#' and \code{lengthx = sapply(x, length)} if \code{x} is a vector, a matrix or a list. 
#' See the examples.
#' 
#' @return  
#' Vectors \code{kmoments} and \code{xmoments} have the following structure 
#' (with a third letter \code{x} added to \code{xmoments}):
#' 
#' \item{ku}{Kurtosis.}
#' \item{ke}{Excess of kurtosis.}
#' \item{sk}{Skewness.}
#' \item{sd}{Standard deviation. Square root of the variance \code{u2}}
#' \item{m1}{Mean.}
#' \item{m2}{Non-central moment of second order.}
#' \item{m3}{Non-central moment of third order.}
#' \item{m4}{Non-central moment of fourth order.}
#' \item{u1}{Central moment of first order. Should be 0.}
#' \item{u2}{Central moment of second order. Variance}
#' \item{u3}{Central moment of third order.}
#' \item{u4}{Central moment of fourth order.}
#' \item{k1}{Cumulant of first order. Should be 0.}
#' \item{k2}{Cumulant of second order.}
#' \item{k3}{Cumulant of third order.}
#' \item{k4}{Cumulant of fourth order.}
#' \item{lh}{Length of x, from which NA values were removed.}
#' \item{......}{.}
#'  
#'   
#' @seealso    
#' \code{\link{pk2pk}}, \code{\link{paramkienerX}}, \code{\link{regkienerLX}}.
#' 
#' @examples
#' 
#' ## Example 1
#' kcmoment(2, c(-1, 1, 6, 9), model = "K2")
#' kcmoment(2, c(-1, 1, 7.2, -0.2/7.2), model = "K3")
#' kcmoment(2, c(-1, 1, 7.2, -0.2), model = "K4")
#' kcmoment(2, c(-1, 1, 6, 7.2, 9, -0.2/7.2, -0.2))
#' kvariance(c(-1, 1, 6, 9))
#' kmoments(c(-1, 1, 6, 9), dgts = 3)
#' 
#' ## Example 2: "K2" and "K7" are preferred input formats for kmoments
#' ## Moments fall at expected parameter values (=> NA).
#' ## apply and direct calculation (= transpose)
#' (mat4 <- matrix(c(rep(0,4), rep(1,4), c(1.9,2.1,3.9,4.1), rep(5,4)),
#'                 nrow = 4, byrow = TRUE, 
#'                 dimnames = list(c("m","g","a","w"), paste0("b",1:4))))
#' round(mat7 <- apply(mat4, 2, pk2pk), 2)
#' round(rbind(mat7, apply(mat7, 2, kmoments)[2:5,]), 2) 
#' round(cbind(t(mat7), kmoments(t(mat7), dgts = 2)[,2:5]), 2) 
#' 
#' ## Example 3: Matrix, timeSeries, xts, zoo + apply 
#' matret    <- 100*diff(log((EuStockMarkets)))
#' (matcoefk <- apply(matret, 2, paramkienerX5, dgts = 2))
#' (matmomk  <- apply(matcoefk, 2, kmoments, lengthx = nrow(matret), dgts = 2))
#' (matmomx  <- apply(matret, 2, xmoments, dgts = 2))
#' rbind(matcoefk, matmomk[2:5,], matmomx[2:5,])
#' 
#' ## Example 4: List + direct calculation = transpose
#' DS   <- getDSdata() ; dimdim(DS) ; class(DS)
#' (pDS <- paramkienerX5(DS, dimnames = FALSE))
#' (kDS <- kmoments(pDS, lengthx = sapply(DS, length), dgts = 3))
#' (xDS <- xmoments( DS, dgts = 3))
#' cbind(pDS, kDS[,2:5], xDS[,2:5])
#' 
#' 
#' @name kmoments
NULL
#' @export
#' @rdname kmoments
kmoments <- function(coefk, model = "K2", lengthx = NA, dgts = NULL, dimnames = FALSE) {

lengthx1 <- lengthx
if (dimdim(lengthx)[1] == 1 && dimdim(lengthx)[2] != 1) {
	if (dimdim(lengthx)[2] != dimdim(coefk)[2]) {
		stop("lengthx is not equal to nrow(coefk)")}
	lengthx1 <- NA
}
z <- switch(dimdimc(coefk), 
	 "1" = .hkmoments(coefk, model=model, lengthx=lengthx1, dgts=dgts),
	 "2" = t(apply(coefk, 1, .hkmoments, model=model, lengthx=lengthx1, dgts=dgts)),
	 "3" = aperm(apply(coefk, c(1,2), .hkmoments, model=model, lengthx=lengthx1, dgts=dgts), c(2,1,3)),
	"-1" = t(sapply(coefk, .hkmoments, model=model, lengthx=lengthx1, dgts=dgts)),
	stop("kmoments cannot handle this format")
	# "numeric" "matrix" "array" "list" "error"
	# "array" params x dates x stocks # aperm dates x params x stocks
	)
if (dimdim(lengthx)[1] == 1 && dimdim(lengthx)[2] != 1) {
	z[,"lh"] <- lengthx
}
if (dimnames) {
	if (dimdim1(z) == 2) {
		dimnames(z) <- list( "STOCKS" = dimnames(z)[[1]],
							"MOMENTS" = dimnames(z)[[2]])
	}
	if (dimdim1(z) == 3) {
		dimnames(z) <- list(  "DATES" = dimnames(z)[[1]],
							"MOMENTS" = dimnames(z)[[2]],
							 "STOCKS" = dimnames(z)[[3]]) 
	}
}
return(z)
}
##
.hkmoments <- function(coefk, model = "K2", lengthx = NA, dgts = NULL) {

coeff  <- pk2pk(coefk, model = model)
z  <- c(
	"ku" = kkurtosis(coeff) , 
	"ke" = kekurtosis(coeff) , 
	"sk" = kskewness(coeff) ,
	"sd" = sqrt(kcmoment(2, coeff)) ,
	"m1" = kmoment(1, coeff) ,
	"m2" = kmoment(2, coeff) ,
	"m3" = kmoment(3, coeff) ,
	"m4" = kmoment(4, coeff) ,
	"u1" = kcmoment(1, coeff) ,
	"u2" = kcmoment(2, coeff) ,
	"u3" = kcmoment(3, coeff) ,
	"u4" = kcmoment(4, coeff) ,
	"k1" = kcmoment(1, coeff) ,
	"k2" = kcmoment(2, coeff) ,
	"k3" = kcmoment(3, coeff) ,
	"k4" = kcmoment(4, coeff) - 3*(kcmoment(2, coeff))^2 ,
	"lh" = lengthx )
names(z) <- c(
	"ku", 
	"ke", 
	"sk",
	"sd",
	"m1",
	"m2",
	"m3",
	"m4",
	"u1",
	"u2",
	"u3",
	"u4",
	"k1",
	"k2",
	"k3",
	"k4",
	"lh" )
z <- if (is.null(dgts)) {z} else {round(z, dgts)}
return(z)
}
#' @export
#' @rdname kmoments
xmoments <- function(x, dgts = NULL, dimnames = FALSE) { 
z <- switch(dimdimc(x), 
	 "1" = .hxmoments(x, dgts=dgts),
	 "2" = t(apply(x, 2, .hxmoments, dgts=dgts)), 
	 "3" = aperm(apply(x, c(1,2), .hxmoments, dgts=dgts), c(2,1,3)), 
	"-1" = t(sapply(x, .hxmoments, dgts=dgts)),
	stop("xmoments cannot handle this format")
	# "numeric" "matrix" "array" "list" "error"
	# "array" params x dates x stocks # aperm dates x params x stocks
	)
if(dimnames) {
	if (dimdim1(z) == 2) {
		dimnames(z) <- list("STOCKS" = dimnames(z)[[1]],
						   "MOMENTS" = dimnames(z)[[2]])
	}
	if (dimdim1(z) == 3) {
		dimnames(z) <- list( "DATES" = dimnames(z)[[1]],
						   "MOMENTS" = dimnames(z)[[2]],
							"STOCKS" = dimnames(z)[[3]]) 
	}
}
return(z)
}
##
.hxmoments <- function(x, dgts = NULL) {
x  <- as.numeric(x[!is.na(x)])
y  <- x - mean(x)
s  <- sd(x)
z  <- c( 
	"kux" = mean(y^4) / s^4 , 
	"kex" =(mean(y^4) / s^4) - 3 , 
	"skx" = mean(y^3) / s^3 ,
	"sdx" = s ,
	"m1x" = mean(x) ,
	"m2x" = mean(x^2) ,
	"m3x" = mean(x^3) ,
	"m4x" = mean(x^4) ,
	"u1x" = mean(y) ,
	"u2x" = mean(y^2) ,
	"u3x" = mean(y^3) ,
	"u4x" = mean(y^4) ,
	"k1x" = mean(y) ,
	"k2x" = mean(y^2) ,
	"k3x" = mean(y^3) ,
	"k4x" = mean(y^4) - 3*(mean(y^2))^2 ,
	"lh"  = length(x) 
	)
z[is.nan(z)] <- NA
z <- if (is.null(dgts)) {z} else {round(z, dgts)}
return(z)
}
#' @export
#' @rdname kmoments
kmoment  <- function(n, coefk, model = "K2", dgts = NULL) {

if (!is.element(n, 1:99)) {stop("n must be an integer 0 < n < 99")}
coeff <- pk2pk(coefk, model = model)
m  <- coeff[1]
g  <- coeff[2]
a  <- coeff[3]
k  <- coeff[4]
w  <- coeff[5] 
if (is.na(min(a, w)) || n >= min(a, w) ) { momk <- NA } else {
	momk <- 0
	for (i in 0:n) {
		for (j in 0:i) { 
			momk <- (momk + choose(n, i) *choose(i, j) 
					 *m^(n-i) *g^(i) *k^(i) *(-1)^(j) 
					 *beta(1 -j/a +(i-j)/w,1 +j/a -(i-j)/w) )  
		} 
	} 
}
names(momk) <- paste0("m", n)
momk <- if (is.null(dgts)) {momk} else {round(momk, dgts)}
return(momk)
}
#' @export 
#' @rdname kmoments
kcmoment  <- function(n, coefk, model = "K2", dgts = NULL) {

if (!is.element(n, 1:99)) {stop("n must be an integer 0 < n < 99")}
coeff <- pk2pk(coefk, model = model)
m  <- coeff[1]
g  <- coeff[2]
a  <- coeff[3]
k  <- coeff[4]
w  <- coeff[5]
if (is.na(min(a, w)) || n >= min(a, w) ) { redcmomk <- NA } else {
	nu <- -beta(1 -1/a, 1 +1/a) + beta(1 -1/w, 1 +1/w)
	redcmomk  <- 0
	for (i in 0:n) {
		for (j in 0:i) { 
		redcmomk <- (redcmomk + choose(n, i) 
		             *choose(i, j) *(-nu)^(n-i) *(-1)^(j) 
					 *beta(1 -j/a +(i-j)/w, 1 +j/a -(i-j)/w) )
		} 
	}
}  
cmomk  <- redcmomk *g^n *k^n
names(cmomk) <- paste0("u", n)
cmomk  <- if (is.null(dgts)) {cmomk} else {round(cmomk, dgts)}
return(cmomk)
}
#' @export 
#' @rdname kmoments
kmean <- function(coefk, model = "K2", dgts = NULL) {
	m1 <- kmoment(1, coefk, model)
	m1 <- if (is.null(dgts)) {m1} else {round(m1, dgts)}
return(m1)
}
#' @export 
#' @rdname kmoments
kstandev <- function(coefk, model = "K2", dgts = NULL) {
	stdev <- sqrt(kcmoment(2, coefk, model))
	names(stdev) <- "sd"
	stdev <- if (is.null(dgts)) {stdev} else {round(stdev, dgts)}
return(stdev)
}
#' @export 
#' @rdname kmoments
kvariance <- function(coefk, model = "K2", dgts = NULL) {
	variance <- kcmoment(2, coefk, model)
	names(variance) <- "u2"
	variance <- if (is.null(dgts)) {variance} else {round(variance, dgts)}
return(variance)
}
#' @export 
#' @rdname kmoments
kskewness <- function(coefk, model = "K2", dgts = NULL) {
	skewk <- kcmoment(3, coefk, model) / kcmoment(2, coefk, model)^(1.5)
	names(skewk) <- "sk"
	skewk <- if (is.null(dgts)) {skewk} else {round(skewk, dgts)}
return(skewk)
}
#' @export 
#' @rdname kmoments
kkurtosis <- function(coefk, model = "K2", dgts = NULL) {
	kurtk <- kcmoment(4, coefk, model) / kcmoment(2, coefk, model)^(2)
	names(kurtk) <- "ku"
	kurtk <- if (is.null(dgts)) {kurtk} else {round(kurtk, dgts)}
return(kurtk)
}
#' @export 
#' @rdname kmoments
kekurtosis <- function(coefk, model = "K2", dgts = NULL) {
	kurtke <- ( kcmoment(4, coefk, model) / kcmoment(2, coefk, model)^(2) ) - 3
	names(kurtke) <- "ke"
	kurtke <- if (is.null(dgts)) {kurtke} else {round(kurtke, dgts)}
return(kurtke)
}


