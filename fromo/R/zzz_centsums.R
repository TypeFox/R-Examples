# This file is part of fromo.
#
# fromo is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fromo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with fromo.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2016.03.30
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav <steven@corecast.io>
# Comments: Steven E. Pav
# Copyright 2016-2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# univariate input#FOLDUP


#' @title centsums Class.
#'
#' @description 
#'
#' An S4 class to store (centered) sums of data, and to support operations on 
#' the same.
#'
#' @details
#'
#' A \code{centsums} object contains a vector value of the data count,
#' the mean, and the \eqn{k}th centered sum, for \eqn{k} up to some
#' maximum order.
#'
#' @slot sums a numeric vector of the sums.
#' @slot order the maximum order.
#'
#' @return An object of class \code{centsums}.
#' @keywords moments
#'
#' @examples 
#' obj <- new("centsums",sums=c(1000,1.234,0.235),order=2)
#'
#' @template etc
#' @template ref-romo
#' @name centsums-class
#' @rdname centsums-class
#' @exportClass centsums
#' @export
setClass("centsums", 
				 representation(sums="numeric",order="numeric"),
				 prototype(sums=c(0.0,0.0,0.0),
									 order=2),
				 validity=function(object) {
					 # ... 
					 # http://www.cyclismo.org/tutorial/R/s4Classes.html
					 if ((!is.null(object@order)) && (length(object@sums) != (1+object@order))) { return("bad dimensionality or order given.") }
					 return(TRUE)
				 }
)
# constructor method documentation
#  
#' @param .Object a \code{centsums} object, or proto-object.
#' @rdname centsums-class
#' @aliases initialize,centsums-class
setMethod('initialize',
					signature('centsums'),
					function(.Object,sums,order=NA_real_) {
						if (is.null(order)) {
							order <- length(sums) - 1
						}
					 	.Object@sums <- sums
					 	.Object@order <- order

						.Object
					})

#'
#' @param sums a numeric vector.
#' @param order the order, defaulting to \code{length(sums)+1}.
#' @name centsums
#' @rdname centsums-class
#' @export
centsums <- function(sums,order=NULL) {
	if (is.null(order)) {
		order <- length(sums) + 1
	}
	retv <- new("centsums", sums=sums, order=order)
	invisible(retv)
}

#' @title Coerce to a centsums object.
#'
#' @description 
#'
#' Convert data to a \code{centsums} object.
#'
#' @details
#'
#' Computes the raw sums on data, and stuffs the results into a 
#' \code{centsums} object.
#'
#' @usage
#'
#' as.centsums(x, order=3, na.rm=TRUE)
#'
#' @param x a numeric, array, or matrix.
#' @param na.rm whether to remove \code{NA}.
#' @inheritParams centsums
#' @return A centsums object.
#' @template etc
#' @examples 
#' set.seed(123)
#' x <- rnorm(1000)
#' cs <- as.centsums(x, order=5)
#' @rdname as.centsums
#' @export as.centsums
as.centsums <- function(x, order=3, na.rm=TRUE) {
	UseMethod("as.centsums", x)
}
#' @rdname as.centsums
#' @export
#' @method as.centsums default
#' @aliases as.centsums
as.centsums.default <- function(x, order=3, na.rm=TRUE) {
	sums <- cent_sums(x, max_order=order, na_rm=na.rm)
	invisible(centsums(sums,order=order))
}

#' @title Accessor methods.
#'
#' @description
#'
#' Access slot data from a \code{centsums} object.
#'
#' @param x a \code{centsums} object.
#' @param type the type of moment to compute.
#' @template etc
#' @name accessor
#' @rdname accessor-methods
#' @aliases sums
#' @exportMethod sums
setGeneric('sums', signature="x", function(x) standardGeneric('sums'))
#' @rdname accessor-methods
#' @aliases sums,centsums-method
setMethod('sums', 'centsums', function(x) x@sums )

# used below
.csums2moments <- function(c_sums,type=c('central','raw','standardized')) {
		# add used_df
		type <- match.arg(type)
		cmoments <- c(c_sums[1],c_sums[2:length(c_sums)] / c_sums[1])

		switch(type,
			raw={
				retv <- cent2raw(cmoments)
			},
			central={ 
				retv <- cmoments[2:length(cmoments)]
				retv[1] <- 0
			},
			standardized={
				retv <- cmoments[2:length(cmoments)]
				retv[1] <- 0.0
				if (length(retv) > 1) {
					if (length(retv) > 2) {
						sigma2 <- retv[2]
						retv[3:length(retv)] <- retv[3:length(retv)] / (sigma2 ^ ((3:length(retv))/2.0))
					}
					retv[2] <- 1.0
				}
			})
			retv
}

#' @rdname accessor-methods
#' @aliases moments
#' @exportMethod moments
setGeneric('moments', function(x,type=c('central','raw','standardized')) standardGeneric('moments'))
#' @rdname accessor-methods
#' @aliases moments,centsums-method
setMethod('moments', signature(x='centsums'),
	function(x,type=c('central','raw','standardized')) {
		# add used_df
		type <- match.arg(type)
		retv <- .csums2moments(x@sums,type)
	})


#' @title concatenate centsums objects.
#' @description 
#'
#' Concatenate centsums objects.
#'
#' @param ... \code{centsums} objects
#' @rdname centsums-concat
#' @seealso join_cent_sums
#' @method c centsums
#' @export
#' @usage \\method{c}{centsums}(...)
c.centsums <- function(...) { 
	.join2 <- function(x,y) {
		x@sums <- join_cent_sums(x@sums,y@sums)
		x
	}
	x <- Reduce(.join2,list(...))
} 
#' @title unconcatenate centsums objects.
#' @description 
#'
#' Unconcatenate centsums objects.
#'
#' @param x a \code{centsums} objects
#' @param y a \code{centsums} objects
#' @seealso unjoin_cent_sums
#' @rdname centsums-unconcat
#' @exportMethod %-%
setGeneric('%-%', function(x,y) standardGeneric('%-%'))
#' @rdname centsums-unconcat
#' @aliases %-%,centsums,centsums-method
setMethod('%-%', signature(x='centsums',y='centsums'),
	function(x,y) {
		x@sums <- unjoin_cent_sums(x@sums,y@sums)
		return(x)
	})

# show#FOLDUP
# 2FIX: add documentation and export
#' @title Show a centsums object.
#'
#' @description 
#'
#' Displays the centsums object.
#'
#' @usage
#'
#' show(object)
#'
#' @param object a \code{centsums} object.
#' @examples 
#' set.seed(123)
#' x <- rnorm(1000)
#' obj <- as.centsums(x, order=5)
#' obj
#' @template etc
#' @name show
#' @rdname show-methods
#' @exportMethod show
#' @aliases show
NULL
#' @rdname show-methods
#' @aliases show,centsums-method
setMethod('show', signature('centsums'), 
# 2FIX: add cumulants?
					function(object) {
						cat('          class:', class(object), '\n')
						cat('    raw moments:', .csums2moments(object@sums,'raw'), '\n')
						cat('central moments:', .csums2moments(object@sums,'central'), '\n')
						cat('    std moments:', .csums2moments(object@sums,'standardized'), '\n')
					})
#UNFOLD
#UNFOLD

# multivariate input#FOLDUP


#' @title centcosums Class.
#'
#' @description 
#'
#' An S4 class to store (centered) cosums of data, and to support operations on 
#' the same.
#'
#' @details
#'
#' A \code{centcosums} object contains a multidimensional array (now only
#' 2-diemnsional), as output by \code{cent_cosums}.
#'
#' @seealso cent_cosums
#' @slot cosums a multidimensional array of the cosums.
#' @slot order the maximum order. ignored for now.
#'
#' @return An object of class \code{centcosums}.
#' @keywords moments
#'
#' @examples 
#' obj <- new("centcosums",cosums=cent_cosums(matrix(rnorm(100*3),ncol=3),max_order=2),order=2)
#'
#' @template etc
#' @template ref-romo
#' @name centcosums-class
#' @rdname centcosums-class
#' @exportClass centcosums
#' @export
setClass("centcosums", 
				 representation(cosums="array",order="numeric"),
				 prototype(cosums=matrix(0,nrow=2,ncol=2),
									 order=2),
				 validity=function(object) {
					 # ... 
					 # http://www.cyclismo.org/tutorial/R/s4Classes.html
					 if ((!is.null(object@order)) && (length(dim(object@cosums)) != object@order)) { return("bad dimensionality or order given.") }
					 if (nrow(object@cosums) != ncol(object@cosums)) { return("must give square cosums for now.") }
					 return(TRUE)
				 }
)
# constructor method documentation
#  
#' @param .Object a \code{centcosums} object, or proto-object.
#' @rdname centcosums-class
#' @aliases initialize,centcosums-class
setMethod('initialize',
					signature('centcosums'),
					function(.Object,cosums,order=NA_real_) {
						if (is.null(order)) {
							order <- 2
						}
					 	.Object@cosums <- cosums
					 	.Object@order <- order

						.Object
					})

#' @param cosums the output of \code{\link{cent_cosums}}, say.
#' @param order the order, defaulting to \code{2}.
#' @name centcosums
#' @rdname centcosums-class
#' @export
centcosums <- function(cosums,order=NULL) {
	if (is.null(order)) {
		order <- length(dim(cosums))
	}
	retv <- new("centcosums", cosums=cosums, order=order)
	invisible(retv)
}

#' @title Coerce to a centcosums object.
#'
#' @description 
#'
#' Convert data to a \code{centcosums} object.
#'
#' @details
#'
#' Computes the raw cosums on data, and stuffs the results into a 
#' \code{centcosums} object.
#'
#' @usage
#'
#' as.centcosums(x, order=2, na.omit=TRUE)
#'
#' @param x a matrix.
#' @param na.omit whether to remove rows with \code{NA}.
#' @inheritParams centcosums
#' @return A centcosums object.
#' @template etc
#' @examples 
#' set.seed(123)
#' x <- matrix(rnorm(100*3),ncol=3)
#' cs <- as.centcosums(x, order=2)
#' @rdname as.centcosums
#' @export as.centcosums
as.centcosums <- function(x, order=2, na.omit=TRUE) {
	UseMethod("as.centcosums", x)
}
#' @rdname as.centcosums
#' @export
#' @method as.centcosums default
#' @aliases as.centcosums
as.centcosums.default <- function(x, order=2, na.omit=TRUE) {
	cosums <- cent_cosums(x, max_order=order, na_omit=na.omit)
	invisible(centcosums(cosums,order=order))
}

#' @title Accessor methods.
#'
#' @description
#'
#' Access slot data from a \code{centcosums} object.
#'
#' @param x a \code{centcosums} object.
#' @param type the type of moment to compute.
#' @template etc
#' @name centcosums-accessor
#' @rdname centcosum-accessor-methods
#' @aliases cosums
#' @exportMethod cosums
setGeneric('cosums', signature="x", function(x) standardGeneric('cosums'))
#' @rdname centcosum-accessor-methods
#' @aliases sums,centcosums-method
setMethod('cosums', 'centcosums', function(x) x@cosums )

# used below
.cosums2comoments <- function(c_sums,type=c('central','raw')) {
		# add used_df
		type <- match.arg(type)
		cmoments <- c(c_sums[1],c_sums[2:length(c_sums)] / c_sums[1])

		switch(type,
			raw={
				retv <- c_sums
				retv[1,1] <- 1
				retv[2:(nrow(retv)),2:(nrow(retv))] <- retv[2:(nrow(retv)),2:(nrow(retv))] + tcrossprod(retv[2:(nrow(retv)),1])
			},
			central={ 
				retv <- c_sums
				retv[1,1] <- 1
				retv[2:(nrow(retv)),1] <- 0
				retv[1,2:(nrow(retv))] <- 0
			})
			retv
}

#' @rdname centcosum-accessor-methods
#' @aliases comoments
#' @exportMethod comoments
setGeneric('comoments', function(x,type=c('central','raw')) standardGeneric('comoments'))
#' @rdname centcosum-accessor-methods
#' @aliases comoments,centcosums-method
setMethod('comoments', signature(x='centcosums'),
	function(x,type=c('central','raw')) {
		# add used_df
		type <- match.arg(type)
		retv <- .cosums2comoments(x@cosums,type)
	})


#' @title concatenate centcosums objects.
#' @description 
#'
#' Concatenate centcosums objects.
#'
#' @param ... \code{centcosums} objects
#' @rdname centcosums-concat
#' @seealso join_cent_cosums
#' @method c centcosums
#' @export
#' @usage \\method{c}{centcosums}(...)
c.centcosums <- function(...) { 
	.join2 <- function(x,y) {
		x@cosums <- join_cent_cosums(x@cosums,y@cosums)
		x
	}
	x <- Reduce(.join2,list(...))
} 
#' @title unconcatenate centcosums objects.
#' @description 
#'
#' Unconcatenate centcosums objects.
#'
#' @param x a \code{centcosums} objects
#' @param y a \code{centcosums} objects
#' @seealso unjoin_cent_cosums
#' @rdname centcosums-unconcat
#' @aliases %-%,centcosums,centcosums-method
setMethod('%-%', signature(x='centcosums',y='centcosums'),
	function(x,y) {
		x@cosums <- unjoin_cent_cosums(x@cosums,y@cosums)
		return(x)
	})

# 2FIX: show a centcosums object
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
