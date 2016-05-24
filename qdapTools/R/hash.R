#' Hash/Dictionary Lookup
#' 
#' \code{hash} - Creates a \href{http://datatable.r-forge.r-project.org/}{\pkg{data.table}} 
#' based hash table for quick hash style dictionary lookup.
#' 
#' @param x A two column dataframe.
#' @param terms A vector of terms to undergo a lookup.
#' @param key The hash key to use.
#' @param missing Value to assign to terms not found in the hash table.
#' @return \code{hash} - Creates a "hash table", a two column \pkg{data.table}. 
#' @seealso 
#' \code{\link[data.table]{setDT}},
#' \code{\link[qdapTools]{hash}}
#' @keywords hash, dictionary, lookup
#' @rdname hash
#' @export
#' @importFrom data.table setkey setDT 
#' @examples
#' ##===================##
#' ## data.table Hashes ##
#' ##===================##
#' (DF <- aggregate(mpg~as.character(carb), mtcars, mean))
#' x <- sample(DF[, 1], 20, TRUE)
#' new.hash <- hash(DF) 
#' x2 <- c(9, 12, x)
#' hash_look(x, new.hash)
#' 
#' x %hl% new.hash
#' x2 %hl% new.hash
#' x2 %hl+% new.hash
#' 
#' ## Create generic functions
#' hfun <- function(x, ...) {
#'     hsh <- hash(x, ...)
#'     function(x, ...) hash_look(x, hsh, ...)
#' }
#' 
#' m <- hfun(DF)
#' m(x)
#' 
#' ##====================##
#' ## Environment Hashes ##
#' ##====================##
#' new.hash2 <- hash_e(DF)
#' 
#' x %hl% new.hash2
#' x2 %hl% new.hash2
#' x2 %hl+% new.hash2
hash <- 
function(x) {

    key <- data.frame(x=x[, 1], y=x[,2], 
        stringsAsFactors = FALSE)

    if (is.factor(key[, 2])) {
        key[, 2] <- as.character(key[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(key[, 2])))
    }

    setDT(key)
    setkey(key, x)        

    class(key) <- c("qdap_hash", "key", class(key))
    attributes(key)[["mode"]] <- FUN
    key
                                                           
}



#' Hash/Dictionary Lookup
#' 
#' \code{hash_look} - Works with a hash table such as is returned from 
#' \code{hash}, to lookup values.
#' 
#' @export
#' @rdname hash
hash_look <- function(terms, key, missing = NA) {

	if (!is.environment(key)) {
		if (is.null(attributes(key)[["mode"]])) {
			return(hash_lookup_helper(terms, key, missing))
		}		
        attributes(key)[["mode"]](hash_lookup_helper(terms, key, missing))
	} else {
		if (is.null(attributes(key)[["mode"]])) {
			return(hash_look_e(terms, key, missing))
		}
		attributes(key)[["mode"]](hash_look_e(terms, key, missing))
	} 
}

#' Hash/Dictionary Lookup
#' 
#' \code{\%hl\%} - A binary operator version of \code{hash_look}.
#'
#' @export
#' @rdname hash
`%hl%` <- function(terms, key) hash_look(terms = terms, key = key)

#' Hash/Dictionary Lookup
#' 
#' \code{\%hl+\%} - A binary operator version of \code{hash_look} 
#' for when \code{missing} is assumed to be \code{NULL}.
#'
#' @export
#' @rdname hash
`%hl+%` <- function(terms, key) hash_look(terms = terms, key = key, missing=NULL)

#' Hash/Dictionary Lookup
#' 
#' \code{\%ha\%} - A deprecated binary operator version of \code{hash_look}.  
#' This will be removed in a subsequent version of \pkg{qdapTools}.  Use 
#' \code{\%hl\%} instead.
#'
#' @param terms A vector of terms to undergo a lookup.
#' @param key The hash key to use.
#' @param key.match Takes one of the following: (1) a two column data.frame of a 
#' match key and reassignment column, (2) a named list of vectors (Note: if 
#' data.frame or named list supplied no key reassign needed) or (3) a single 
#' vector match key.
#' @export
#' @rdname Deprecated
`%ha%` <- function(terms, key) {
    .Deprecated(msg = paste("`%ha%` is deprecated.  Please use `%hl%` instead."), 
    	old = as.character(sys.call(sys.parent()))[1L])
	
	hash_look(terms = terms, key = key)
}

### Helper function
#' @importFrom data.table setDT
hash_lookup_helper <- function(terms, key, missing = NA) {
	
    terms <- data.frame(x=terms, stringsAsFactors = FALSE)
    setDT(terms)
 
    out <- key[terms][[2]]

    if (!is.null(missing) && is.na(missing)) return(out)
    if (!is.null(missing) && !is.na(missing)) {
        hits <- which(is.na(out))
        out[hits] <- missing
        return(out)
    }

    if (is.null(missing)) {
        hits <- which(is.na(out))
        out[hits] <- terms[[1]][hits]
        return(out)
    }

}
