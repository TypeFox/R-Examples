

#' @include d_logishp.R




#' @title Length and Dimensions of Vector, Matrix, Array, Data.Frame, List
#'
#' @description
#' Length and dimensions of vector, matrix, array, data.frame and list. 
#' A friendly version of \code{dim} that returns the true dimension 
#' rather than the sometimes unexpected \code{NULL} value. 
#' The number of dimensions appears first, then the length in each dimension. 
#' A special case is list: the list length (number of items) is turned into a 
#' negative integer and the length/dimension of each item is either positive 
#' if the item is a vector, matrix, array or data.frame or negative if the item 
#' is itself a list. Only the first level of the list is explored. 
#' \code{dimdim1} and \code{dimdimc} return the first item of \code{dimdim}, 
#' thus the true dimension, either as an integer or as a character 
#' (in this latest case, always \code{"-1"} for lists). 
#' Note: Some problems may occur with S4 objects like 
#' \code{dimdim(qualityTools::fracDesign(k = 3, gen = "C = AB"))}.
#' 
#' @param        x    vector, matrix, array, data.frame, list.
#' 
#' @examples
#' 
#' require(timeSeries)
#' 
#' dimdim(NULL) ; dimdim(NA) ; dimdim(NaN)
#' dimdim(Inf) ; dimdim(TRUE) ; dimdim(FALSE)
#' dimdim(11:39)
#' dimdim(LETTERS[1:8])
#' dimdim(matrix(1:60, ncol=5))
#' dimdim(extractData())
#' dimdim(as.data.frame(extractData()))
#' dimdim(data.frame(X=1:2, Y=1:4, Z=LETTERS[1:8]))
#' dimdim(array(1:240, c(8,6,5)))
#' dimdim(array(1:240, c(4,2,6,5)))
#' dimdim(getDSdata())
#' dimdim(zData)
#' dimdim(xData)
#' dimdim(tData)
#' dimdim1(tData)
#' dimdimc(tData)
#' 
#' @export
#' @name dimdim
dimdim    <- function(x) { 
z  <- 	if (is.list(x) & !is.data.frame(x)) { 
				c(-length(x), sapply(x, dimdim1))
		} else { 
			if (!is.null(dim(x))) { 
				c(length(dim(x)), dim(x)) 
			} else { 
				if (max(NCOL(x), ncol(x), na.rm=TRUE) == 1) { 
					c(NCOL(x), NROW(x)) 
				} else { 
                    c(2, max(NROW(x), nrow(x), na.rm=TRUE), 
					     max(NCOL(x), ncol(x), na.rm=TRUE)) 
				}
			}
		}
names(z) <- NULL
return(z)
}
#' @export
#' @rdname dimdim
dimdim1 <- function(x) { dimdim(x)[1] }
#' @export
#' @rdname dimdim
dimdimc <- function(x) {
	ddim  <- dimdim(x)[1] 
	z     <- if (ddim < 0) {"-1"} else {as.character(ddim)}
return(z)
}


#' @title Local Conversion Functions Between Kiener Distribution Parameters
#'
#' @description
#' Conversion functions between parameters \code{a}, \code{k}, \code{w}, 
#' \code{d}, \code{e} used in Kiener distributions K2, K3 and K4.
#' 
#' @param        a    a numeric value.
#' @param        k    a numeric value.
#' @param        w    a numeric value.
#' @param        d    a numeric value.
#' @param        e    a numeric value.
#' 
#' @details
#' \code{a} (alpha) is the left tail parameter, 
#' \code{w} (omega) is the right tail parameter, 
#' \code{d} (delta) is the distortion parameter, 
#' \code{e} (epsilon) is the eccentricity parameter. 
#' \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} and 
#' describes a global tail parameter. 
#' They are defined by: 
#' \deqn{ aw2k(a, w) = k = 2 / (1/a + 1/w) = \frac{2}{\frac{1}{a} +\frac{1}{w}}  }
#' \deqn{ aw2d(a, w) = d = (-1/a + 1/w) / 2 = \frac{-\frac{1}{a} +\frac{1}{w}}{2} }
#' \deqn{ aw2e(a, w) = e = (a - w) / (a + w) = \frac{a-w}{a+w} }
#' \deqn{ kd2a(k, d) = a = 1 / ( 1/k - d) = \frac{1}{\frac{1}{k} - d} }
#' \deqn{ kd2w(k, d) = w = 1 / ( 1/k + d) = \frac{1}{\frac{1}{k} + d} }
#' \deqn{ ke2a(k, e) = a = k / (1 - e) = \frac{k}{1-e} }
#' \deqn{ ke2w(k, e) = w = k / (1 + e) = \frac{k}{1+e} }
#' \deqn{ ke2d(k, e) = d = e / k = \frac{e}{k} }
#' \deqn{ kd2e(k, d) = e = k * d }
#' \deqn{ de2k(k, e) = k = e / d = \frac{e}{d} }
#'
#' @seealso 
#' The asymmetric Kiener distributions K2, K3, K4:  
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}
#'
#' @examples
#' 
#' aw2k(4, 6); aw2d(4, 6); aw2e(4, 6)
#' outer(1:6, 1:6, aw2k)
#' 
#' @export
#' @name aw2k
           aw2k <- function(a, w) { 2 / (1/a + 1/w) }
#' @export
#' @rdname aw2k
           aw2d <- function(a, w) { (-1/a + 1/w) / 2 }
#' @export
#' @rdname aw2k
           aw2e <- function(a, w) { (a - w) / (a + w) }
#' @export
#' @rdname aw2k
           ad2e <- function(a, d) { d*a / ( 1 + d*a) }
#' @export
#' @rdname aw2k
           ad2k <- function(a, d) { a / ( 1 + d*a) }
#' @export
#' @rdname aw2k
           ad2w <- function(a, d) { 1 / ( 2*d + 1/a) }
#' @export
#' @rdname aw2k
           ae2d <- function(a, e) { e / (1 - e) / a }
#' @export
#' @rdname aw2k
           ae2k <- function(a, e) { (1 - e) * a }
#' @export
#' @rdname aw2k
           ae2w <- function(a, e) { a * (1 - e) / (1 + e) }
#' @export
#' @rdname aw2k
           ak2d <- function(a, k) { (a - k) / a / k }
#' @export
#' @rdname aw2k
           ak2e <- function(a, k) { (a - k) / a }
#' @export
#' @rdname aw2k
           ak2w <- function(a, k) { 1 / ( 2/k - 1/a) }
#' @export
#' @rdname aw2k
           de2a <- function(d, e) { e / d / (1 - e) }
#' @export
#' @rdname aw2k
           de2k <- function(d, e) { e / d }
#' @export
#' @rdname aw2k
           de2w <- function(d, e) { e / d / (1 + e) }
#' @export
#' @rdname aw2k
           dk2a <- function(d, k) { 1 / (1/k - d) }
#' @export
#' @rdname aw2k
           dk2e <- function(d, k) { d * k }
#' @export
#' @rdname aw2k
           dk2w <- function(d, k) { 1 / (1/k + d) }
#' @export
#' @rdname aw2k
           dw2a <- function(d, w) { -1 / ( 2*d - 1/w) }
#' @export
#' @rdname aw2k
           dw2e <- function(d, w) { d*w / ( 1 - d*w) }
#' @export
#' @rdname aw2k
           dw2k <- function(d, w) { w / ( 1 - d*w) }
#' @export
#' @rdname aw2k
           ek2a <- function(e, k) { k / (1 - e) }
#' @export
#' @rdname aw2k
           ek2d <- function(e, k) { e / k }
#' @export
#' @rdname aw2k
           ek2w <- function(e, k) { k / (1 + e) }
#' @export
#' @rdname aw2k
           ew2a <- function(e, w) { w * (1 + e) / (1 - e) }
#' @export
#' @rdname aw2k
           ew2d <- function(e, w) { e / (1 + e) / w }
#' @export
#' @rdname aw2k
           ew2k <- function(e, w) { (1 + e) * w }
#' @export
#' @rdname aw2k
           kd2a <- function(k, d) { 1 / ( 1/k - d) }
#' @export
#' @rdname aw2k
           kd2e <- function(k, d) { k * d }
#' @export
#' @rdname aw2k
           kd2w <- function(k, d) { 1 / ( 1/k + d) }
#' @export
#' @rdname aw2k
           ke2a <- function(k, e) { k / (1 - e) }
#' @export
#' @rdname aw2k
           ke2d <- function(k, e) { e / k }
#' @export
#' @rdname aw2k
           ke2w <- function(k, e) { k / (1 + e) }
#' @export
#' @rdname aw2k
           kw2a <- function(k, w) { 1 / ( 2/k - 1/w) }
#' @export
#' @rdname aw2k
           kw2d <- function(k, w) { (k - w) / w / k }
#' @export
#' @rdname aw2k
           kw2e <- function(k, w) { (k - w) / w }



#' @title Global Conversion Function Between Kiener Distribution Parameters
#'
#' @description
#' A conversion function between Kiener distribution parameters 
#' \code{K1(m, g, k)}, \code{K2(m, g, a, w)},
#' \code{K3(m, g, k, d)} and \code{K4(m, g, k, e)} to and from
#' \code{coefk = c(m, g, a, k, w, d, e)} extracted from \code{\link{regkienerLX}} 
#' and \code{\link{paramkienerX}}. 
#'  
#' @param     coefk    vectors of numeric of length 3, 4 or 7.
#' @param     model    character. Either "K1", "K2", "K3", "K4", "K7".
#' @param     to       character. Either "K1", "K2", "K3", "K4", "K7".
#' @param     dgts     integer. The rounding applied to the output.	
#' 
#' @details
#' Kiener distributions use the following parameters, some of them being redundant. 
#' See also \code{\link{aw2k}} for the formulas and 
#' the conversion between parameters:
#' \itemize{
#'   \item{ \code{m} (mu) is the median of the distribution,. }
#'   \item{ \code{g} (gamma) is the scale parameter. }
#'   \item{ \code{a} (alpha) is the left tail parameter. } 
#'   \item{ \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} 
#'          and describes a global tail parameter. }
#'   \item{ \code{w} (omega) is the right tail parameter. } 
#'   \item{ \code{d} (delta) is the distortion parameter. }
#'   \item{ \code{e} (epsilon) is the eccentricity parameter. }
#' }
#' 
#' \code{pk2pk()} performs the conversion between the various representation, from and to:
#' \itemize{
#'   \item{ "K1" : \code{kiener1(m, g, k)}    }
#'   \item{ "K2" : \code{kiener2(m, g, a, w)} }
#'   \item{ "K3" : \code{kiener3(m, g, k, d)} }
#'   \item{ "K4" : \code{kiener4(m, g, k, e)} }
#'   \item{ "K7" : \code{c(m, g, a, k, w, d, e)} }
#' }
#' 
#' \code{coefk} can take any of the above form. When length(coefk) is 4, 
#' \code{model = "K2", "K3" or "K4"} is required to differentiate the three models.
#' When length(coefk) is 3 or 7, recognition is automatic and  
#' \code{model = "K1" or "K7"} is ignored. The vector is assumed to be correct 
#' and there is no check of the consistency between the 
#' parameters \code{a, k, w, d} and \code{e}.
#' 
#' The output may be any of the above forms. Default is \code{"K7" = c(m, g, a, k, w, d, e)} 
#' which is \code{coefk} provided by the regression function \code{\link{regkienerLX}}
#' or the parameter estimation function \code{\link{paramkienerX}}. It is widely in many plots.
#' 
#' An integer rounding parameter is provided trough \code{dgts}. Default is no rounding.
#'
#' @seealso 
#' Local conversion functions \code{\link{aw2k}}, 
#' Kiener distributions K1, K2, K3 and K4: \code{\link{kiener1}}, 
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}
#'
#' @examples
#'
#' ## Example 1
#' c2 <- c(1, 2, 3, 5)
#' pk2pk(c2, model = "K2", to = "K1") # loose the asymmetry.
#' pk2pk(c2, model = "K2", to = "K2")
#' pk2pk(c2, model = "K2", to = "K3")
#' pk2pk(c2, model = "K2", to = "K4")
#' pk2pk(c2, model = "K2", to = "K4")
#' (c7 <- pk2pk(c2, model = "K2", to = "K7", dgts = 3))
#' pk2pk(c7, model = "K7", to = "K2")
#' 
#' ## Example 2 ("K2" to "K7")
#' (mat4 <- matrix( c(rep(0,9), rep(1,9), seq(0.5,4.5,0.5), seq(1,5,0.5)), 
#'          nrow = 4, byrow = TRUE, dimnames = list(c("m","g","a","w"), paste0("b",1:9))))
#' (mat7 <- round(apply(mat4, 2, pk2pk), 3))
#' 
#' 
#' @export
#' @name pk2pk	   
pk2pk <- function(coefk, model = "K2", to = "K7", dgts = NULL) {
model    <- toupper(model) # for compatibility with old versions
to       <- toupper(to)
coeff    <- if (is.matrix(coefk)) { 
	switch( as.character(ncol(coefk)) , 
	"3"  = cbind(coefk[,1], coefk[,2], coefk[,3], coefk[,3], coefk[,3], 
	             rep(0, nrow(coefk)), rep(0, nrow(coefk))),
	"4"  = switch( model , 
			"K2"  = cbind(coefk[,1], coefk[,2], coefk[,3], aw2k(coefk[,3], coefk[,4]), 
					  coefk[,4], aw2d(coefk[,3], coefk[,4]), aw2e(coefk[,3], coefk[,4])) ,
			"K3"  = cbind(coefk[,1], coefk[,2], kd2a(coefk[,3], coefk[,4]), coefk[,3],  
					  kd2w(coefk[,3], coefk[,4]), coefk[,4], kd2e(coefk[,3], coefk[,4])) ,
			"K4"  = cbind(coefk[,1], coefk[,2], ke2a(coefk[,3], coefk[,4]), coefk[,3],  
					  ke2w(coefk[,3], coefk[,4]), ke2d(coefk[,3], coefk[,4]), coefk[,4]) , 
			stop("when ncol(coefk) = 4, model must be either K2, K3 or K4.")
			) ,
	"7"  = coefk ,
	stop("coefk is of wrong size. ncol(coefk) must be a matrix of 3, 4 or 7 columns.")
	)} else { 
	switch( as.character(length(coefk)) , 
	"3"  = c(coefk[1], coefk[2], coefk[3], coefk[3], coefk[3], 
	         0, 0) , 
	"4"  = switch( model, 
			"K2"  = c(coefk[1], coefk[2], coefk[3], aw2k(coefk[3], coefk[4]), 
					  coefk[4], aw2d(coefk[3], coefk[4]), aw2e(coefk[3], coefk[4])) , 
			"K3"  = c(coefk[1], coefk[2], kd2a(coefk[3], coefk[4]), coefk[3],  
					  kd2w(coefk[3], coefk[4]), coefk[4], kd2e(coefk[3], coefk[4])) ,
			"K4"  = c(coefk[1], coefk[2], ke2a(coefk[3], coefk[4]), coefk[3],  
					  ke2w(coefk[3], coefk[4]), ke2d(coefk[3], coefk[4]), coefk[4]) , 
			stop("when length(coefk) = 4, model must be either K2, K3 or K4.")
			) ,
	"7"  = coefk ,
	stop("coefk is of wrong size. length(coefk) must be 3, 4 or 7.")
	)}
zcoefk <- if (is.matrix(coeff)) { 
		switch( to , 
			"K1"  = coeff[,c(1,2,4)], 
			"K2"  = coeff[,c(1,2,3,5)], 
			"K3"  = coeff[,c(1,2,4,6)],
			"K4"  = coeff[,c(1, 2, 4, 7)], 
			"K7"  = coeff, 
			stop("to must be either K1, K2, K3, K4 or K7")
			)} else {
		switch( to , 
			"K1"  = coeff[c(1,2,4)], 
			"K2"  = coeff[c(1,2,3,5)], 
			"K3"  = coeff[c(1,2,4,6)],
			"K4"  = coeff[c(1, 2, 4, 7)], 
			"K7"  = coeff, 
			stop("to must be either K1, K2, K3, K4 or K7")			
			)}
if (is.matrix(zcoefk)) { 
		colnames(zcoefk) <- switch( to , 
			"K1"  = c("m", "g", "k") , 
			"K2"  = c("m", "g", "a", "w") , 
			"K3"  = c("m", "g", "k", "d") ,
			"K4"  = c("m", "g", "k", "e") , 
			"K7"  = c("m", "g", "a", "k", "w", "d", "e") 
			)} else {
		names(zcoefk) <- switch( to , 
			"K1"  = c("m", "g", "k") , 
			"K2"  = c("m", "g", "a", "w") , 
			"K3"  = c("m", "g", "k", "d") ,
			"K4"  = c("m", "g", "k", "e") , 
			"K7"  = c("m", "g", "a", "k", "w", "d", "e") 
			)}
if (!is.null(dgts)) {zcoefk <- round(zcoefk, dgts)}
return(zcoefk)
}




#' @title Round Coefk
#'
#' @description
#' Round coefk parameters in a standard manner or in a special manner, 
#' the latest being useful to display nice matrix or data.frame.
#' 
#' @param    coefk      numeric, matrix or data.frame representing
#'                      parameters \code{c(m,g,a,k,w,d,e)}.
#' @param    dgts       integer. The number of rounded digits. 
#' @param    parnames   boolean. Output displayed with or without parameter names.
#' 
#' @details
#' For \code{dgts} between 1 and 9, rounding is done in the standard way
#' and all parameters are rounded at the same number of digits.
#' 
#' For \code{dgts} between 10 and 27, rounding of parameters 
#' \code{c(m,g,a,k,w,d,e)} is done in the following way: 
#' \itemize{
#'   \item{ dgts = 10 : c(0, 0, 1, 1, 1, 3, 2)}
#'   \item{ dgts = 11 : c(1, 1, 1, 1, 1, 3, 2)}
#'   \item{ dgts = 12 : c(2, 2, 1, 1, 1, 3, 2)}
#'   \item{ dgts = 13 : c(3, 3, 1, 1, 1, 3, 2)}
#'   \item{ dgts = 14 : c(4, 4, 1, 1, 1, 3, 2)}
#'   \item{ dgts = 15 : c(5, 5, 1, 1, 1, 3, 2)}
#'   \item{ dgts = 16 : c(0, 0, 2, 2, 2, 3, 2)}
#'   \item{ dgts = 17 : c(1, 1, 2, 2, 2, 3, 2)}
#'   \item{ dgts = 18 : c(2, 2, 2, 2, 2, 3, 2)}
#'   \item{ dgts = 19 : c(3, 3, 2, 2, 2, 3, 2)}
#'   \item{ dgts = 20 : c(4, 4, 2, 2, 2, 3, 2)}
#'   \item{ dgts = 21 : c(5, 5, 2, 2, 2, 3, 2)}
#'   \item{ dgts = 22 : c(0, 0, 3, 3, 3, 4, 3)}
#'   \item{ dgts = 23 : c(1, 1, 3, 3, 3, 4, 3)}
#'   \item{ dgts = 24 : c(2, 2, 3, 3, 3, 4, 3)}
#'   \item{ dgts = 25 : c(3, 3, 3, 3, 3, 4, 3)}
#'   \item{ dgts = 26 : c(4, 4, 3, 3, 3, 4, 3)}
#'   \item{ dgts = 27 : c(5, 5, 3, 3, 3, 4, 3)}
#' }
#' 
#' @examples     
#' 
#' mat   <- matrix(runif(35), ncol=7) 
#' coefk <- mat[1,]
#' 
#' roundcoefk(coefk, dgts = 2, parnames = FALSE)
#' roundcoefk(coefk, dgts = 15)
#' roundcoefk(mat, dgts = 15)
#' 
#' 
#' @export
#' @name roundcoefk
roundcoefk <- function(coefk, dgts = NULL, parnames = TRUE) {
	z <- if (is.matrix(coefk) || is.data.frame(coefk)) {
				t(apply(coefk, 1, .hroundcoefk, dgts = dgts, parnames = parnames))
			} else { 
				.hroundcoefk(coefk, dgts = dgts, parnames = parnames)
			}
return(z)
}

##
.hroundcoefk <- function(coefk, dgts, parnames) {
	z <- if (is.null(dgts)) { 
			coefk 
		} else { 
		   switch(as.character(dgts),
				"10" = round(coefk, digits = c(0, 0, 1, 1, 1, 3, 2)),
				"11" = round(coefk, digits = c(1, 1, 1, 1, 1, 3, 2)),
				"12" = round(coefk, digits = c(2, 2, 1, 1, 1, 3, 2)),
				"13" = round(coefk, digits = c(3, 3, 1, 1, 1, 3, 2)),
				"14" = round(coefk, digits = c(4, 4, 1, 1, 1, 3, 2)),
				"15" = round(coefk, digits = c(5, 5, 1, 1, 1, 3, 2)),
				"16" = round(coefk, digits = c(0, 0, 2, 2, 2, 3, 2)),
				"17" = round(coefk, digits = c(1, 1, 2, 2, 2, 3, 2)),
				"18" = round(coefk, digits = c(2, 2, 2, 2, 2, 3, 2)),
				"19" = round(coefk, digits = c(3, 3, 2, 2, 2, 3, 2)),
				"20" = round(coefk, digits = c(4, 4, 2, 2, 2, 3, 2)),
				"21" = round(coefk, digits = c(5, 5, 2, 2, 2, 3, 2)),
				"22" = round(coefk, digits = c(0, 0, 3, 3, 3, 4, 3)),
				"23" = round(coefk, digits = c(1, 1, 3, 3, 3, 4, 3)),
				"24" = round(coefk, digits = c(2, 2, 3, 3, 3, 4, 3)),
				"25" = round(coefk, digits = c(3, 3, 3, 3, 3, 4, 3)),
				"26" = round(coefk, digits = c(4, 4, 3, 3, 3, 4, 3)),
				"27" = round(coefk, digits = c(5, 5, 3, 3, 3, 4, 3)),
				round(coefk, digits = dgts) ) 
		}
	names(z) <- if (parnames) { c("m","g","a","k","w","d","e")} else { NULL }
return(z)
}


#' @title Generate a list of vectors of characters from a vector of probabilities
#' 
#' @description
#' Generate vector of characters from a vector of probabilities, replacing 
#' \code{0.} by letters:
#' \itemize{
#'   \item{ \code{p.} : probability. }
#'   \item{ \code{q.} : quantile. }
#'   \item{ \code{VaR.} : Value-at-Risk, positive in most cases. } 
#'   \item{ \code{c.} : corrective tail coefficient = (q - m) / (q_logistic_function - m). }
#'   \item{ \code{ltm.} : left tail mean (signed ES on the left tail, usually negative). } 
#'   \item{ \code{rtm.} : right tail mean (signed ES on the right tail, usually positive). }
#'   \item{ \code{dtmq.} : (p<=0.5 left, p>0.5 right) tail mean minus quantile. }
#'   \item{ \code{ES.} : expected shortfall, positive in most cases. }
#'   \item{ \code{h.} : corrective ES  = (ES - m) / (ES_logistic_function - m). }
#'   \item{ \code{desv.} : ES - VaR, usually positive. } 
#'   \item{ \code{l.} : quantile of the tangent logistic function. }
#'   \item{ \code{dl.} : quantile - quantile_logistic_function. }
#'   \item{ \code{g.} : quantile of the Laplace-Gauss function. } 
#'   \item{ \code{dg.} : quantile - quantile_Laplace_Gauss_function. }
#' }
#' 
#' 
#' , \code{q.}, \code{VaR.}, \code{c.}, \code{ltm.},
#' \code{rtm.}, \code{ES.}, \code{h.}, \code{l.}, \code{dl.}, \code{g.}, \code{dg.}.
#' The result is a list of vectors.
#' 
#' @param    probak    a vector of ordered probabilities with 0 and 1 excluded.
#'
#' @seealso 
#' Probabilities: \code{\link{pprobs0}}
#' 
#' @examples     
#' getnamesk(pprobs1)
#' getnamesk(pprobs8)
#' 
#' @export
#' @name getnamesk
getnamesk <- function(probak = pprobs2) {

if ((sum(probak <= 0) + sum(probak >= 1)) > 0) { stop("probak is not within (0, 1).") }
if (!checkquantiles(probak)) { stop("probak is not ordered.") }
nprobak  <- character(length(probak))
for (i in 1:length(probak)) { 
	nprobak[i] <- sub("0.", "p.", format(probak[i], nsmall=2, scientific=FALSE)) 
	} 
namesk <- list(
	  probak = probak,
	 nprobak = nprobak,
	 nquantk = sub("p.", "q.",   nprobak),
	   nvark = sub("p.", "VaR.", nprobak),
	 nctailk = sub("p.", "c.",   nprobak),
	   nltmk = sub("p.", "ltm.", nprobak),
	   nrtmk = sub("p.", "rtm.", nprobak),
	  ndtmqk = sub("p.", "dtmq.",nprobak),
	    nesk = sub("p.", "ES.",  nprobak),
	   nhesk = sub("p.", "h.",   nprobak),
	  ndesvk = sub("p.", "desv.",nprobak),
	 nlogisk = sub("p.", "l.",   nprobak),
	ndlogisk = sub("p.", "dl.",  nprobak),
	 ngaussk = sub("p.", "g.",   nprobak),
	ndgaussk = sub("p.", "dg.",  nprobak)
	)
return(namesk)
}
