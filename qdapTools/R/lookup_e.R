#' Hash Table/Dictionary Lookup
#'  
#' \code{lookup_e} - Environment based hash table useful for large vector lookups.
#' 
#' @param terms A vector of terms to undergo a lookup_e.
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
#' \code{\link[base]{new.env}}, \code{\link[qdapTools]{lookup}},
#' @keywords dictionary, hash, lookup
#' @export
#' @rdname lookup_e
#' @examples
#' lookup_e(1:5, data.frame(1:4, 11:14))
#' 
#' ## Retain original values for missing
#' lookup_e(1:5, data.frame(1:4, 11:14), missing=NULL)
#' 
#' lookup_e(LETTERS[1:5], data.frame(LETTERS[1:5], 100:104))
#' lookup_e(LETTERS[1:5], factor(LETTERS[1:5]), 100:104)
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
#' lookup_e(1:10, codes)
#' 
#' ## Supply a single vector to key.match and key.reassign
#' 
#' lookup_e(mtcars$carb, sort(unique(mtcars$carb)),
#'     c("one", "two", "three", "four", "six", "eight"))
#' 
#' lookup_e(mtcars$carb, sort(unique(mtcars$carb)),
#'     seq(10, 60, by=10))
#' 
#' ## %le%, a binary operator version of lookup
#' 1:5 %le% data.frame(1:4, 11:14)
#' 1:10 %le% codes
#' 
#' 1:12 %le% codes
#' 1:12 %le+% codes
lookup_e <-
function(terms, key.match, key.reassign=NULL, missing = NA) {

    UseMethod("lookup_e", key.match)

}

#' @export
#' @method lookup_e matrix
#' @rdname lookup_e
lookup_e.matrix <-
function(terms, key.match, key.reassign=NULL, missing = NA) {

    key.match <- data.frame(key.match, stringsAsFactors = FALSE)

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
    }
    mode.out <- mode(key.match[, 2])    
    key.match[, 1] <- as.character(key.match[, 1])

    hits <- which(!is.na(match(terms, key.match[, 1])))
    x <- rep(ifelse(is.null(missing), NA, missing), length(terms))
	
    KEY <- hash_e(key.match, mode.out = mode.out)   
    x[hits] <- recoder(terms[hits], envr = KEY)

    if (is.null(missing)) { 
    	keeps <- which(is.na(x))
        x[keeps] <- terms[keeps]
        x
    }   
    x
}

#' @export
#' @method lookup_e data.frame
#' @rdname lookup_e
lookup_e.data.frame <-
function(terms, key.match, key.reassign=NULL, missing = NA) {

    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
    }
    mode.out <- mode(key.match[, 2])    
    key.match[, 1] <- as.character(key.match[, 1])

    hits <- which(!is.na(match(terms, key.match[, 1])))
    x <- rep(ifelse(is.null(missing), NA, missing), length(terms))
	
    KEY <- hash_e(key.match, mode.out = mode.out)   
    x[hits] <- recoder(terms[hits], envr = KEY)

    if (is.null(missing)) { 
    	keeps <- which(is.na(x))
        x[keeps] <- terms[keeps]
        x
    }   
	x
}

#' @export
#' @method lookup_e list
#' @rdname lookup_e
lookup_e.list <-
function(terms, key.match, key.reassign=NULL, missing = NA) {

    key.match <- list2df(key.match) 
    if (is.factor(key.match[, 2])) {
        key.match[, 2] <- as.character(key.match[, 2])
    }
    mode.out <- mode(key.match[, 2])    
    key.match[, 1] <- as.character(key.match[, 1])

    hits <- which(!is.na(match(terms, key.match[, 1])))
    x <- rep(ifelse(is.null(missing), NA, missing), length(terms))
	
    KEY <- hash_e(key.match, mode.out = mode.out)   
    x[hits] <- recoder(terms[hits], envr = KEY)

    if (is.null(missing)) { 
    	keeps <- which(is.na(x))
        x[keeps] <- terms[keeps]
        x
    }    
    x
}

#' @export
#' @method lookup_e numeric
#' @rdname lookup_e
lookup_e.numeric <-
function(terms, key.match, key.reassign=NULL, missing = NA) {
     
    mode.out <- mode(key.reassign)    
    DF <- data.frame(as.character(key.match), key.reassign, 
        stringsAsFactors = FALSE)   

    hits <- which(!is.na(match(terms, DF[, 1])))
    x <- rep(ifelse(is.null(missing), NA, missing), length(terms))
	
    KEY <- hash_e(DF, mode.out = mode.out)   
    x[hits] <- recoder(terms[hits], envr = KEY)

    if (is.null(missing)) { 
    	keeps <- which(is.na(x))
        x[keeps] <- terms[keeps]
        x
    }   
    x
}


#' @export
#' @method lookup_e factor
#' @rdname lookup_e
lookup_e.factor <-
function(terms, key.match, key.reassign=NULL, missing = NA) {

    key.reassign <- as.character(key.reassign)
        
    mode.out <- mode(key.reassign)    
    DF <- data.frame(as.character(key.match), key.reassign, 
        stringsAsFactors = FALSE)   

    hits <- which(!is.na(match(terms, DF[, 1])))
    x <- rep(ifelse(is.null(missing), NA, missing), length(terms))
	
    KEY <- hash_e(DF, mode.out = mode.out)   
    x[hits] <- recoder(terms[hits], envr = KEY)

    if (is.null(missing)) { 
    	keeps <- which(is.na(x))
        x[keeps] <- terms[keeps]
        x
    }   
    x
}

#' @export
#' @method lookup_e character
#' @rdname lookup_e
lookup_e.character <-
function(terms, key.match, key.reassign=NULL, missing = NA) {

    mode.out <- mode(key.reassign)    	
    DF <- data.frame(as.character(key.match), key.reassign, 
        stringsAsFactors = FALSE)   

    hits <- which(!is.na(match(terms, DF[, 1])))
    x <- rep(ifelse(is.null(missing), NA, missing), length(terms))
	
    KEY <- hash_e(DF, mode.out = mode.out)   
    x[hits] <- recoder(terms[hits], envr = KEY)

    if (is.null(missing)) { 
    	keeps <- which(is.na(x))
        x[keeps] <- terms[keeps]
        x
    }   
    x
}


#' Hash/Dictionary Lookup
#' 
#' \code{\%le\%} - A binary operator version of \code{lookup_e} 
#' for when \code{key.match} is a data.frame or named list.
#'
#' @export
#' @rdname lookup_e
`%le%` <- function(terms, key.match) lookup_e(terms = terms, key.match = key.match)

#' Hash/Dictionary Lookup
#' 
#' \code{\%le+\%} - A binary operator version of \code{lookup_e} 
#' for when \code{key.match} is a data.frame or named list and \code{missing} is
#' assumed to be \code{NULL}.
#'
#' @export
#' @rdname lookup_e
`%le+%` <- function(terms, key.match) {
    lookup_e(terms = terms, key.match = key.match, missing = NULL)
}


