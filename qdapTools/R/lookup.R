#' Hash Table/Dictionary Lookup
#'  
#' \code{lookup} - \href{http://datatable.r-forge.r-project.org/}{\pkg{data.table}} 
#' based hash table useful for large vector lookups.
#' 
#' @param terms A vector of terms to undergo a lookup.
#' @param key.match Takes one of the following: (1) a two column data.frame of a 
#' match key and reassignment column, (2) a named list of vectors (Note: if 
#' data.frame or named list supplied no key reassign needed) or (3) a single 
#' vector match key.
#' @param key.reassign A single reassignment vector supplied if key.match is 
#' not a two column data.frame/named list.
#' @param missing Value to assign to terms not matching the key.match.  If set 
#' to \code{NULL} the original values in \code{terms} corresponding to the 
#' missing elements are retained.
#' @return Outputs A new vector with reassigned values.
#' @seealso 
#' \code{\link[data.table]{setDT}},
#' \code{\link[qdapTools]{hash}}
#' @keywords dictionary, hash, lookup
#' @export
#' @rdname lookup
#' @examples
#' ## Supply a dataframe to key.match
#' 
#' lookup(1:5, data.frame(1:4, 11:14))
#' 
#' ## Retain original values for missing 
#' lookup(1:5, data.frame(1:4, 11:14), missing=NULL) 
#' 
#' lookup(LETTERS[1:5], data.frame(LETTERS[1:5], 100:104))
#' lookup(LETTERS[1:5], factor(LETTERS[1:5]), 100:104)
#' 
#' ## Supply a named list of vectors to key.match
#' 
#' codes <- list(
#'     A = c(1, 2, 4), 
#'     B = c(3, 5),
#'     C = 7,
#'     D = c(6, 8:10)
#' )
#' 
#' lookup(1:10, codes)
#' 
#' ## Supply a single vector to key.match and key.reassign
#' 
#' lookup(mtcars$carb, sort(unique(mtcars$carb)),        
#'     c("one", "two", "three", "four", "six", "eight")) 
#'     
#' lookup(mtcars$carb, sort(unique(mtcars$carb)),        
#'     seq(10, 60, by=10))
#'   
#' ## %l%, a binary operator version of lookup
#' 1:5 %l% data.frame(1:4, 11:14)
#' 1:10 %l% codes
#' 
#' 1:12 %l% codes
#' 1:12 %l+% codes
#'   
#' (key <- data.frame(a=1:3, b=factor(paste0("l", 1:3))))
#' 1:3 %l% key
#' 
#' ##Larger Examples
#' key <- data.frame(x=1:2, y=c("A", "B"))
#' big.vec <- sample(1:2, 3000000, TRUE)
#' out <- lookup(big.vec, key)
#' out[1:20]
#' 
#' ## A big string to recode with variation
#' ## means a bigger dictionary
#' recode_me <- sample(1:(length(LETTERS)*10), 10000000, TRUE)
#' 
#' ## Time it
#' tic <- Sys.time()  
#' 
#' output <- recode_me %l% split(1:(length(LETTERS)*10), LETTERS)
#' difftime(Sys.time(), tic)
#' 
#' ## view it
#' sample(output, 100)
lookup <-
function(terms, key.match, key.reassign=NULL, missing = NA) {

    UseMethod("lookup", key.match)

}


#' @export
#' @method lookup list
#' @rdname lookup
lookup.list <-
function (terms, key.match, key.reassign = NULL, missing = NA) {
    key.match <- list2df(key.match, "x", "y")

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(key.match[, 2])))
    }
	
	output <- lookup_helper(terms, key.match, missing)

	if(attributes(output)[["missing"]]) return(FUN(output))

    out_warn <- tryCatch({
        FUN(output)
    }, warning = function(w) {
        TRUE
    }, finally = {
        FALSE
    })
	
    if(length(out_warn) == 1 && !isTRUE(out_warn)) return(FUN(output))
	
    attributes(output) <- NULL
    output
}

#' @export
#' @method lookup data.frame
#' @rdname lookup
lookup.data.frame <-
function (terms, key.match, key.reassign = NULL, missing = NA) {
	
    key.match <- data.frame(x=key.match[, 1], y=key.match[,2], 
        stringsAsFactors = FALSE)

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(key.match[, 2])))
    }

    output <- lookup_helper(terms, key.match, missing)

    if(attributes(output)[["missing"]]) return(FUN(output))

    out_warn <- tryCatch({
        FUN(output)
    }, warning = function(w) {
        TRUE
    }, finally = {
        FALSE
    })
	
    if(length(out_warn) == 1 && !isTRUE(out_warn)) return(FUN(output))
	
    attributes(output) <- NULL
    output

}

#' @export
#' @method lookup matrix 
#' @rdname lookup
lookup.matrix <-
function (terms, key.match, key.reassign = NULL, missing = NA) {
    key.match <- data.frame(x=key.match[, 1], y=key.match[,2], 
        stringsAsFactors = FALSE)

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(key.match[, 2])))
    }

	output <- lookup_helper(terms, key.match, missing)

	if(attributes(output)[["missing"]]) return(FUN(output))

    out_warn <- tryCatch({
        FUN(output)
    }, warning = function(w) {
        TRUE
    }, finally = {
        FALSE
    })
	
    if(length(out_warn) == 1 && !isTRUE(out_warn)) return(FUN(output))
	
    attributes(output) <- NULL
    output

}

#' @export
#' @method lookup numeric 
#' @rdname lookup
lookup.numeric <-
function(terms, key.match, key.reassign, missing = NA) {
     
    key.match <- data.frame(x=key.match, y=key.reassign, 
        stringsAsFactors = FALSE)

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(key.match[, 2])))
    }

	output <- lookup_helper(terms, key.match, missing)

	if(attributes(output)[["missing"]]) return(FUN(output))

    out_warn <- tryCatch({
        FUN(output)
    }, warning = function(w) {
        TRUE
    }, finally = {
        FALSE
    })
	
    if(length(out_warn) == 1 && !isTRUE(out_warn)) return(FUN(output))
	
    attributes(output) <- NULL
    output
}

#' @export
#' @method lookup factor
#' @rdname lookup
lookup.factor <-
function(terms, key.match, key.reassign, missing = NA) {
     
    key.match <- data.frame(x=key.match, y=key.reassign, 
        stringsAsFactors = FALSE)

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(key.match[, 2])))
    }

	output <- lookup_helper(terms, key.match, missing)

	if(attributes(output)[["missing"]]) return(FUN(output))

    out_warn <- tryCatch({
        FUN(output)
    }, warning = function(w) {
        TRUE
    }, finally = {
        FALSE
    })
	
    if(length(out_warn) == 1 && !isTRUE(out_warn)) return(FUN(output))
	
    attributes(output) <- NULL
    output
}

#' @export
#' @method lookup character 
#' @rdname lookup
lookup.character <-
function(terms, key.match, key.reassign, missing = NA) {
     
    key.match <- data.frame(x=key.match, y=key.reassign, 
        stringsAsFactors = FALSE)

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(key.match[, 2])))
    }

	output <- lookup_helper(terms, key.match, missing)

	if(attributes(output)[["missing"]]) return(FUN(output))

    out_warn <- tryCatch({
        FUN(output)
    }, warning = function(w) {
        TRUE
    }, finally = {
        FALSE
    })
	
    if(length(out_warn) == 1 && !isTRUE(out_warn)) return(FUN(output))
	
    attributes(output) <- NULL
    output
}


#' @importFrom data.table setkey setDT
lookup_helper <- function(terms, key, missing = NA) {

	x <- i.y <- NULL
	
    terms <- data.frame(x=terms, stringsAsFactors = FALSE)
    key <- data.table::data.table(key[c("x", "y")])
    setDT(terms)
 
    setkey(key, x)
    out <- key[terms][[2]]
    attributes(out) <- list(missing = TRUE)
	
    if (!is.null(missing) && is.na(missing)) return(out)
    if (!is.null(missing) && !is.na(missing)) {
        hits <- which(is.na(out))
        out[hits] <- missing
        return(out)
    }

    if (is.null(missing)) {
        hits <- which(is.na(out))
        out[hits] <- terms[[1]][hits]
        attributes(out) <- list(missing = FALSE)
        return(out)
    }

}



#' Hash/Dictionary Lookup
#' 
#' \code{\%l\%} - A binary operator version of \code{lookup} 
#' for when \code{key.match} is a data.frame or named list.
#'
#' @export
#' @rdname lookup
`%l%` <- function(terms, key.match) lookup(terms = terms, key.match = key.match)

#' Hash/Dictionary Lookup
#' 
#' \code{\%l+\%} - A binary operator version of \code{lookup} 
#' for when \code{key.match} is a data.frame or named list and \code{missing} is
#' assumed to be \code{NULL}.
#'
#' @export
#' @rdname lookup
`%l+%` <- function(terms, key.match) {
    lookup(terms = terms, key.match = key.match, missing = NULL)
}

#' Hash/Dictionary Lookup
#' 
#' \code{\%lc\%} - A binary operator version of \code{lookup} 
#' for when \code{key.match} is a data.frame or named list and all arguments are 
#' converted to character.
#'
#' @export
#' @rdname lookup
`%lc%` <- function(terms, key.match) {
	
	terms <- as.character(terms)
	if (!is.data.frame(key.match) && is.list(key.match)) {
	    key.match <- list2df(key.match, "x", "y")	
	}
	key.match[] <- lapply(key.match, as.character)
	
    lookup(terms = terms, key.match = key.match)
}

#' Hash/Dictionary Lookup
#' 
#' \code{\%lc+\%} - A binary operator version of \code{lookup} 
#' for when \code{key.match} is a data.frame or named list, \code{missing} is
#' assumed to be \code{NULL}, and all arguments are converted to character.
#'
#' @export
#' @rdname lookup
`%lc+%` <- function(terms, key.match) {
	
	terms <- as.character(terms)
	if (!is.data.frame(key.match) && is.list(key.match)) {
	    key.match <- list2df(key.match, "x", "y")	
	}
	key.match[] <- lapply(key.match, as.character)
	
    lookup(terms = terms, key.match = key.match, missing = NULL)
}





