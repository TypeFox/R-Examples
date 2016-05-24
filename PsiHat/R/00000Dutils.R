
#-----------------------
is_vide<-is_nothing<-function(object){length(object)==0}
is_any <- Is <- function(object, class2){any(sapply(class2, is, object = object))}
assert.is <- function(object, class2, text = ""){
	stopifnot(is.character(class2))
	if(missing(object))
		stop(paste(class2, "object missing in assert.is", text))
	stopifnot(length(class2) >= 1)

	if(!is_any(object = object, class2 = class2))
	{
		warning(paste("got ", class(object), sep = ""))
		message("got ", class(object), " when one of these classes was required:")
		print(class2)
		stop(text)
	}
}
assert.are <- function(object, class2, ...){
	assert.is(object, "list")
	if(!is_nothing(object)){
		for(obj in object)
			assert.is(object = obj, class2 = class2, ...)
	}
}
are <- function(object, class2){all(sapply(object, is_any, class2 = class2))}
is_err <- function(object){is(object, "try-error")}

isInteger <- function(x){floor(x) == x}


#-----------------------------------------
	

setGeneric("is_prob", function(object,...) standardGeneric("is_prob"))
setMethod("is_prob",signature(object = "numeric" ), function(object,tolerance = 1e-3, ...){
         boo <- all(object >= 0 - tolerance & object <= 1 + tolerance, ...)
	 !is.na(boo) && boo})
are_prob <- function(object, ...){all(sapply(object, is_prob, ...))}

#--------------------
default <- function(object, name, verbose, return.value = object)
{
	if(missing(verbose))
		verbose <- TRUE
	if(verbose)
	{
		if(missing(name))
			name <- "actual argument"
		if(is.function(name))
			name <- name(object)
		stopifnot(length(name) == 1 && is.character(name))
		prefix <- paste(name, "was set to default value of ")
		suffix <- paste(" on ", date(), ".", sep = "")
		object.name <- if(isS4(object) && any(c("annotation", "ann") %in% slotNames(object)))
		{
			object.nam <- try(paste("object of annotation '", object@annotation, "'", sep = ""))#nannotation(object)
			if(!is(object.nam, "character"))
			{ message("bad object.nam from object0"); browser()}
			object.nam
		}
		else
			object
		if(length(object.name) == 1 && (is(object.name, "character") || is(object.name, "numeric")))
			message(prefix, object.name, suffix)
		else
		{
			cat("\n", prefix)
			print(object.name)
			cat(paste(suffix, "\n\n"))
		}
	}
	return.value
}



printInvalid <- function(object, ...)
{
	mess <- paste("Invalid ", class(object), ..., " on ", date(), ".", sep = "")
	message(mess)
	warning(mess)
}

printGeneric <- function(object, ...)
{
	message("Generic ", if(is.character(object)) object else class(object), ..., " on ", date(), ".")
}



#------------------------
setClass("Numeric", representation("numeric"))
setValidity("Numeric", function(object)
{
	num.ok <- all(is.na(object)) || all(object >= 0, na.rm = TRUE)
	if(!num.ok)
	{ printInvalid(object); browser()}
	num.ok
})
setAs(from = "Numeric", to = "numeric", function(from)
{
  x <- from@.Data
  names(x) <- names(from)
  x
})

#
setClass("scalar", representation("numeric"))
setValidity("scalar", function(object)
{
	length(object) %in% c(0, 1) # changed 16 April 2008
})

setClass("Scalar", representation("scalar"))
setValidity("Scalar", function(object)
{
	length(object) == 0 || is.na(object) || object >= 0 # changed 16 April 2008
})
setIs("Scalar", "Numeric")
#
setClass("Matrix", representation("matrix"))
setValidity("Matrix", function(object)
{
	all(as.numeric(object) >= 0, na.rm = TRUE)
})
##
setClassUnion("Vector", c("numeric", "character", "logical"))
#
setGeneric("nNumeric", function(object,...) standardGeneric("nNumeric"))#Numeric <- function(object, ...){message('Numeric generic'); browser()}
setMethod("nNumeric", signature(object = "numeric"), function(object)
{
  x <- new("Numeric", as.numeric(object))
  names(x) <- names(x@.Data) <- names(object)
  x
})

setGeneric("nscalar", function(object,...) standardGeneric("nscalar"))#scalar <- function(object, ...){message("scalar generic"); browser()}
setMethod("nscalar", signature(object = "scalar"), function(object)
{
	as(object, "scalar")
})
setMethod("nscalar", signature(object = "ANY"), function(object)
{
	nscalar(as.numeric(object))
})
setMethod("nscalar", signature(object = "numeric"), function(object) # changed 16 April 2008
{
#  if(length(object) == 0)
#    stop("cannot convert vector of length 0 to scalar")
#  else
  if(length(object) > 1)
  {
#    message("extra vector elements will be lost in coersion of object to scalar")
    print(stats(object))
    stop("extra vector elements prevents coersion of object to scalar")
  	object <- object[1]
  }
  new("scalar", object)
})







setGeneric("nScalar", function(object,...) standardGeneric("nScalar"))#Scalar <- function(object, ...){message("Scalar generic"); browser()}
setMethod("nScalar", signature(object = "ANY"), function(object) # changed 16 April 2008
{
#  if(length(object) == 0)
#    stop("cannot convert vector of length 0 to Scalar")
#  else if(length(object) > 1)
#    warning("extra vector elements lost in coersion to Scalar")
  Sca <- nscalar(object)
  if(length(Sca) == 1 && !is.na(Sca) && Sca < 0)
    stop("cannot convert negative number to Scalar")
  new("Scalar", Sca)
})


setGeneric("nMatrix", function(object,...) standardGeneric("nMatrix"))#
setMethod("nMatrix", signature(object = "matrix"), function(object)
{
  new("Matrix", object)
})
setMethod("nMatrix", signature(object = "numeric"), function(object)
{
  nNumeric(object) # [sic]
})

##
setMethod("[", signature(x = "Numeric", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  nNumeric(as(x, "numeric")[i])
})
setReplaceMethod("[", signature(x = "Numeric", i = "ANY", j = "missing"), function(x, i, j, value) # showMethods("[<-")
{
	vec <- as(x, "numeric")
	vec[i] <- value
#	tr <- try(vec[i] <- as(value, "numeric"))
#	if(is(tr, "try-error") || length(vec) != length(x))
#	{ message("bad replacement of ", class(x)); browser()}
	nNumeric(vec)
})
setMethod("[", signature(x = "Matrix", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  mat <- as(x, "matrix")
  mat <- mat[i, , drop = FALSE]
  nMatrix(mat)
})
setMethod("[", signature(x = "Matrix", i = "missing", j = "ANY"), function(x, i, j, drop)
{
  mat <- as(x, "matrix")
  mat <- mat[, j, drop = FALSE]
  nMatrix(mat)
})
setMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY"), function(x, i, j, drop)
{
  mat <- as(x, "matrix")
  mat <- mat[i, j, drop = FALSE]
  nMatrix(mat)
})
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "missing"), function(x, i, j, value) # showMethods("[<-")
{
#	cat("replacing with "); print(value)
	mat <- as(x, "matrix")
	tr <- try(mat[i, ] <- value)
	if(is(tr, "try-error"))
	{ message("bad replacement"); browser()}
	nMatrix(mat)
})


##-------------------------------
setGeneric("sameNames", function(object,...) standardGeneric("sameNames"))#sameNames <- function(object, ...){message("generic sameNames"); browser()}
setMethod("sameNames", signature(object = "ANY"), function(object, ...)
{
	arglis <- list(...) # c(list(object), list(...))
	if("order.sensitive" %in% names(arglis))
	{
		order.sensitive <- arglis$order.sensitive
		arglis <- arglis[names(arglis) != "order.sensitive"]
		if("order.sensitive" %in% names(arglis))
		{ message("bad arglis in sameNames"); browser()}
	}
	else
		order.sensitive <- TRUE
	na.rm <- FALSE
	verbose <- FALSE
	nam0 <- names(object) # arglis[[1]])
	nam <- try(if(order.sensitive) nam0 else sort(nam0))
	if(is(nam, "try-error"))
	{
		message("nam err")
		browser()
	}
	if(verbose)
		print(nam)
	all(sapply(arglis, function(elem)
	{
	  elem.nam <- names(elem)
	  if(verbose)
	  	print(elem.nam)
	  if(!order.sensitive)
	  	elem.nam <- sort(elem.nam)
		identical(elem.nam, nam)
	}), na.rm = na.rm)
})


sameLengths <- function(...)
{
	lis <- list(...)
	all(length(lis[[1]]) == sapply(lis, length))
}
compatible <- function(...){sameLengths(...) && sameNames(...)}
#
hsm <- function(x, na.rm) # half-sample mode of a vector [D. R. Bickel and R. Frhwirth (contributed equally), "On a Fast, Robust Estimator of the Mode: Comparisons to Other Robust Estimators with Applications," Computational Statistics and Data Analysis 50, 3500-3530 (2006)]
{
	stopifnot(is.numeric(x))
	if(missing(na.rm))
		na.rm <- FALSE
	if(na.rm)
		x <- x[!is.na(x)]
  y <- sort(x);
  while (length(y)>=4)
  {
    m <- ceiling(length(y)/2);
    w.min <- y[length(y)]-y[1];
    for(i in 1:(length(y)-m+1))
    {
      w <- y[i+m-1]-y[i];
      if(w<=w.min)
      {
        w.min <- w;
        j <- i
      }
    }
    if(w==0)
      y <- y[j]
    else
      y <- y[j:(j+m-1)]
  }
  if(length(y) == 3)
  {
    z <- 2*y[2]-y[1]-y[3];
    if(!is.finite(z))
    {
      print('ERROR: z is not finite; x, y, and z follow:');
      print(x);
      print(y);
      print(z);
    }
    if(z < 0)
      mean(y[1:2])
    else if(z > 0)
      mean(y[2:3])
    else
      y[2]
  }
  else
    mean(y)
}


print_stats <- function(object, name, ...){message(name, ":"); print(stats(object, ...))}
stats <- function(object, ...)
{
	if(is(object, "numeric"))
		object <- as(object, "numeric")
	summ <- try(summary(object, ...))
	if(is_err(summ))
		summ <- c(mean = mean(object), median = median(object), n.na = sum(is.na(object)))
	if(is(object, "numeric"))
	{
		nam <- c("n.finite", "size", "mode", names(summ))
		vec <- c(n.finite = sum(is.finite(object)), size = length(object), mode = hsm(object), as.numeric(summ))
		names(vec) <- nam
		vec
	}
	else
		summ

}

