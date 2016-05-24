rbind.zoo <- function(..., deparse.level = 1)
{  

  args <- Filter(Negate(is.null), list(...))
  indexes <- do.call("c", unname(lapply(args, index)))

  my.table <- function(x) {
    ix <- ORDER(x)
    x <- x[ix]
    table(MATCH(x,x))
  }
  if(max(my.table(indexes)) > 1L) stop("indexes overlap")

  if(any(sapply(args, function(x) is.null(dim(x)) && length(x) == 0L && length(index(x)) > 0L)))
    stop("zero-length vectors with non-zero-length index are not allowed")  

  ncols <- sapply(args, NCOL)
  if(!all(ncols == ncols[1L])) stop("number of columns differ")
  
  ## process colnames() if any
  nams <- lapply(args, colnames)
  namsNULL <- sapply(nams, is.null)
  if(all(namsNULL)) namsOK <- TRUE else {
    if(sum(namsNULL) > 0L) namsOK <- FALSE else {
      nam1 <- nams[[1L]]
      namsID <- sapply(nams, function(x) identical(x, nam1))
      if(all(namsID)) namsOK <- TRUE else {
        namsSORT <- sapply(nams, function(x) identical(sort(x), sort(nam1)))
	if(!all(namsSORT)) namsOK <- FALSE else {
	  namsOK <- TRUE
	  for(i in which(!namsID)) args[[i]] <- args[[i]][,nam1]
	}
      }
    }
  }
  if(!namsOK) warning("column names differ")

  ## collect data
  argsdata <- lapply(args, coredata)

  ## (special case: rbinding of vectors with one-column matrices)
  nulldim <- sapply(argsdata, function(a) is.null(dim(a)))
  if(ncols[1L] == 1L) {
    if(nulldim[1] & any(!nulldim)) {
      argsdata <- lapply(argsdata, function(a) if(is.null(dim(a))) a else a[,1, drop = TRUE])
      nulldim <- rep(TRUE, length(nulldim))
    }
    if(!nulldim[1] & any(nulldim)) {
      argsdata <- lapply(argsdata, function(a) if(is.null(dim(a))) as.matrix(a) else a)
      nulldim <- rep(FALSE, length(nulldim))
    }
  }

  if((ncols[1L] > 1L) | !all(nulldim))
    rval <- zoo(do.call("rbind", argsdata), indexes)
  else
    rval <- zoo(do.call("c", argsdata), indexes)

  freq <- if(!("zooreg" %in% unlist(sapply(args, class)))) NULL
            else {
	      freq <- c(frequency(rval), unlist(sapply(args, frequency)))
	      if((length(freq) == (length(args)+1L)) && 
	         identical(all.equal(max(freq)/freq, round(max(freq)/freq)), TRUE))
		 max(freq) else NULL
	    }
  if(!is.null(freq)) {
    attr(rval, "frequency") <- freq
    class(rval) <- c("zooreg", class(rval))
  }
  return(rval)
}

c.zoo <- function(...) {
    rbind.zoo(...)
}

cbind.zoo <- function(..., all = TRUE, fill = NA, suffixes = NULL, drop = FALSE)
{
  merge.zoo(..., all = all, fill = fill, suffixes = suffixes, retclass = "zoo", drop = drop)
}


merge.zoo <- function(..., all = TRUE, fill = NA, suffixes = NULL, check.names = FALSE, retclass = c("zoo", "list", "data.frame"), drop = TRUE)
{
    if (!is.null(retclass)) retclass <- match.arg(retclass)
    # cl are calls to the args and args is a list of the arguments
    cl <- as.list(match.call())
    cl[[1]] <- cl$all <- cl$fill <- cl$retclass <- cl$suffixes <- cl$check.names <- cl$drop <- NULL
    args <- list(...)

    parent <- parent.frame()

    is.plain <- function(x) 
	all(class(x) %in% c("array", "integer", "numeric", "factor", "matrix", "logical"))

    is.scalar <- function(x) is.plain(x) && length(x) == 1

    # ensure all ... plain args are of length 1 or have same NROW as arg 1
    stopifnot(all(sapply(args, function(x) is.zoo(x) || !is.plain(x) ||
      (is.plain(x) && (NROW(x) == NROW(args[[1]]) || is.scalar(x))))))

    scalars <- sapply(args, is.scalar)

    if(!is.zoo(args[[1]])) args[[1]] <- as.zoo(args[[1]])
    for(i in seq_along(args))
        if (is.plain(args[[i]]))  
            args[[i]] <- zoo(args[[i]], index(args[[1]]), attr(args[[1]], "frequency"))
	else if (!is.zoo(args[[i]]))
            args[[i]] <- as.zoo(args[[i]])

    ## retain frequency	if all series have integer multiples of the same frequency
    ## and at least one of the original objects is a "zooreg" object	
    freq <- if(!("zooreg" %in% unlist(sapply(args, class)))) NULL
        else {
	  freq <- unlist(sapply(args, frequency))
	  if((length(freq) == length(args)) && 
	     identical(all.equal(max(freq)/freq, round(max(freq)/freq)), TRUE))
	     max(freq) else NULL
	}

    # use argument names if suffixes not specified
    if (is.null(suffixes)) {
        makeNames <- function(l) {
            nm <- names(l)
            fixup <- if (is.null(nm)) 
                seq_along(l)
            else nm == ""
            dep <- sapply(l[fixup], function(x) deparse(x)[1])
            if (is.null(nm)) 
                return(dep)
            if (any(fixup)) 
                nm[fixup] <- dep
            nm
        }
        suffixes <- makeNames(as.list(substitute(list(...)))[-1])
    }
    if (length(suffixes) != length(cl)) {
        warning("length of suffixes and does not match number of merged objects")
        suffixes <- rep(suffixes, length.out = length(cl))
    }

    # extend all to a length equal to the number of args
    all <- rep(as.logical(all), length.out = length(cl))

    ## check indexes:
    indexlist <- lapply(args, index)
    ## 1. for non-unique entries
    index_duplicates <- function(x) length(unique(MATCH(x, x))) < length(x)
    if(any(sapply(indexlist, index_duplicates)))
      stop("series cannot be merged with non-unique index entries in a series")
    ## 2. for differing classes
    indexclasses <- sapply(indexlist, function(x) class(x)[1])
    if (!all(indexclasses == indexclasses[1])) 
        warning(paste("Index vectors are of different classes:", 
            paste(indexclasses, collapse = " ")))

    # fn to get the unique elements in x, in sorted order, using only
    # [, match, length and order
    sort.unique <- function(x) {
        ix <- MATCH(x, x) == seq_len(length(x))
	x <- x[ix]
        ix <- ORDER(x)
        x[ix]
    }

    # fn to get intersection of each element in list of lists
    intersect.list <- function(list) { 
        my.table <- function(x) {
	   ix <- ORDER(x)
           x <- x[ix]
           table(MATCH(x, x))
	}
	union <- do.call("c", unname(list))
	sort.unique(union)[which(my.table(union) == length(list))]
    }
    indexintersect <- intersect.list(indexlist)

    # get the indexes of the final answer which is the union of
    # all indexes of args corresponding to all=TRUE with the intersection
    # of all indexes
    indexunion <- do.call("c", unname(indexlist[all]))
    
    indexes <-  if(is.null(indexunion)) indexintersect
      else sort.unique(c(indexunion, indexintersect))
    # previously, we used to do:
    # if (is.null(indexunion)) indexunion <- do.call("c", indexlist)[0]
    # indexes <- sort.unique(c(indexunion, indexintersect))

    ## check whether resulting objects still got the same frequency
    if(!is.null(freq)) {
      freq <- c(frequency(zoo(,indexes)), freq)
      freq <- if((length(freq) == 2) && identical(all.equal(max(freq)/freq, round(max(freq)/freq)), TRUE))
        max(freq) else NULL
    }

    # the f function does the real work
    # it takes a zoo object, a, and fills in a matrix corresponding to
    # indexes with the values in a. ret.zoo is TRUE if it is to return
    # a zoo object.  If ret.zoo is FALSE it simply returns with the matrix
    # just calculated.  
    # match0 is convenience wrapper for MATCH with nomatch=0 default
    match0 <- function(a, b, nomatch = 0, ...) MATCH(a, b, nomatch = nomatch, ...)
    f <- if (any(all)) {
       function(a, ret.zoo = TRUE) {
        if (length(a) == 0 && length(dim(a)) == 0)
	   return(if(ret.zoo) {
	            rval <- zoo(, indexes)
	            attr(rval, "frequency") <- freq
	            if(!is.null(freq)) class(rval) <- c("zooreg", class(rval))
		    rval
		  } else numeric())
        z <- matrix(fill, length(indexes), NCOL(a))
	if (length(dim(a)) > 0)
           z[match0(index(a), indexes), ] <- a[match0(indexes, index(a)), , drop = FALSE]        
        else {
           z[match0(index(a), indexes), ] <- a[match0(indexes, index(a))]
           z <- z[, 1, drop=TRUE]
        }
 	if (ret.zoo) {
	  z <- zoo(z, indexes)
	  attr(z, "oclass") <- attr(a, "oclass")
	  attr(z, "levels") <- attr(a, "levels")
	  attr(z, "frequency") <- freq
	  if(!is.null(freq)) class(z) <- c("zooreg", class(z))
	}
	return(z)
      }
    
    } else {
    # if all contains only FALSE elements then the following f is used
    # instead of the prior f for performance purposes.  If all contains
    # only FALSE then the resulting index is the intersection of the
    # index of each argument so we can just return a[index] or a[index,].
    # Also if we are not to return a zoo object then unclass it prior to return.
      function(a, ret.zoo = TRUE) {
	if (!ret.zoo) class(a) <- NULL
	if (length(dim(a)) == 0) {
		if (length(a) == 0) {
		   rval <- if(ret.zoo) zoo(, indexes) else numeric()
		} else
		   rval <- as.zoo(a[match0(indexes, attr(a, "index"))])
	} else
		rval <- as.zoo(a[match0(indexes, attr(a, "index")), , drop=FALSE])
        if(is.zoo(rval) && !is.null(freq)) {
	  attr(rval, "frequency") <- freq
	  class(rval) <- unique(c("zooreg", class(rval)))
	}
	return(rval)
      }
    }

    # if retclass is NULL do not provide a return value but instead
    # update each argument that is a variable, i.e. not an expression,
    # in place.  
    if (is.null(retclass)) {
        for(vn in cl) {
           if (is.name(vn))
           tryCatch(
	     eval(substitute(v <- f(v), list(f = f, v = vn)), parent), 
	     condition = function(x) {}
           )
        }
	invisible(return(NULL))
    } 

    # apply f to each arg, put result of doing this on all args in list rval
    # and then cbind that list together to produce the required matrix
    rval <- lapply(args, f, ret.zoo = retclass %in% c("list", "data.frame"))
    ## have commented this next line out.  Is it needed?
    # for(i in which(scalars)) rval[[i]] <- rval[[i]][] <- zoo(coredata(rval[[i]])[1], index(rval[[1]]), freq)
    names(rval) <- suffixes
    if (retclass == "list") { 
	return(rval)
    }
    if (retclass == "data.frame") {
      ## transform list to data.frame
      ## this is simple if all list elements are vectors, but with
      ## matrices a bit more effort seems to be needed:
      charindex <- index2char(index(rval[[1]]), frequency = freq)
      nam1 <- names(rval)
      rval <- lapply(rval, as.list)
      todf <- function(x) {
        class(x) <- "data.frame"
        attr(x, "row.names") <- charindex
        return(x)
      }
      rval <- lapply(rval, todf)
      ## name processing
      nam2 <- sapply(rval, function(z) 1:NCOL(z))
      for(i in 1:length(nam2)) nam2[[i]] <- paste(names(nam2)[i], nam2[[i]], sep = ".")
      nam1 <- unlist(ifelse(sapply(rval, NCOL) > 1, nam2, nam1))
      rval <- do.call("cbind", rval)
      names(rval) <- nam1
      ## turn zoo factors into plain factors
      is.zoofactor <- function(x) !is.null(attr(x, "oclass")) && attr(x, "oclass") == "factor"
      for(i in 1:NCOL(rval)) if(is.zoofactor(rval[,i])) rval[,i] <- coredata(rval[,i])
      return(rval)
    }
    # remove zero length arguments
    rval <- rval[sapply(rval, function(x) length(x) > 0)]
    # if there is more than one non-zero length argument then cbind them
    # Note that cbind will create matrices, even when given a single vector, 
    # so - if drop=TRUE - we do not to use it in the single vector case.
    rval <- if (length(rval) > 1L | (length(rval) == 1L & !drop))
	do.call("cbind", rval)
    else if (length(rval) > 0L)
	rval[[1]]
    # return if vector since remaining processing is only for column names
    if (length(dim(rval)) == 0L) {
      # fixed bug: coredata was missing
      rval <- zoo(coredata(rval), indexes)
      attr(rval, "frequency") <- freq
      if(!is.null(freq)) class(rval) <- c("zooreg", class(rval))
      return(rval)
    }

    # processing from here on is to compute nice column names
    if (length(unlist(sapply(args, colnames))) > 0) {
        fixcolnames <- function(a) {
            # if (length(a) == 0) 
            #   return(NULL)
            if (length(dim(a)) ==0) {
                if (length(a) == 0) return(NULL) else return("")
            } else {
				if (ncol(a) == 0) return(NULL)
                rval <- colnames(a)
                if (is.null(rval)) {
                  rval <- paste(1:NCOL(a), suffixes[i], sep = ".")
                }
                else {
                  rval[rval == ""] <- as.character(which(rval == ""))
                }
                return(rval)
            }
        }
        zoocolnames <- lapply(args, fixcolnames)
        zcn <- unlist(zoocolnames)
        fixme <- lapply(zoocolnames, function(x) x %in% zcn[duplicated(zcn)])
        f <- function(i) {
            rval <- zoocolnames[[i]]
            rval[rval == ""] <- suffixes[i]
            rval
        }
        zoocolnames <- lapply(seq_along(args), f)
        f <- function(i) ifelse(fixme[[i]], paste(zoocolnames[[i]], 
            suffixes[i], sep = "."), zoocolnames[[i]])
        if (any(duplicated(unlist(zoocolnames)))) 
            zoocolnames <- lapply(seq_along(args), f)
    } else {
        fixcolnames <- function(a) {
            if (length(a) == 0) 
                return(NULL)
            if (NCOL(a) < 2) 
                return("")
            else return(paste(".", 1:NCOL(a), sep = ""))
        }
        zoocolnames <- lapply(args, fixcolnames)
        zoocolnames <- lapply(seq_along(args), function(i) 
		if (!is.null(zoocolnames[[i]])) # NULL returned if false
			paste(suffixes[i], zoocolnames[[i]], sep = ""))
    }
	zoocolnames <- unlist(zoocolnames)
	colnames(rval) <- if (check.names) make.names(make.unique(zoocolnames))
		else if (ncol(rval) == length(zoocolnames)) zoocolnames else 
		colnames(rval)
    # rval <- zoo(rval, indexes)
    rval <- zoo(coredata(rval), indexes)
    attr(rval, "frequency") <- freq
    if(!is.null(freq)) class(rval) <- c("zooreg", class(rval))
    return(rval)
}


