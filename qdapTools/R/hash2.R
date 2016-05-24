#' Hash/Dictionary Lookup
#' 
#' \code{hash_e} - Creates a new environment for quick hash style dictionary lookup.
#' 
#' @param mode.out The type of output (column 2) expected (e.g.,
#' \code{"character"}, \code{"numeric"}, etc.)
#' @return \code{hash_e} - Creates a "hash table", a two column \code{data.frame} 
#' in its own environment.  
#' @author \code{hash_e} - Bryan Goodrich and Tyler Rinker <tyler.rinker@@gmail.com>.
#' @seealso \code{\link[base]{environment}}
#' @references \url{http://www.talkstats.com/showthread.php/22754-Create-a-fast-dictionary}
#' @rdname hash
#' @export
hash_e <- 
function(x, mode.out = "numeric") {
    hash_help_e(x, mode.out = mode.out)                                                                   
}

hash_help_e <- function(x, mode.out) {
	
    if (is.factor(x[, 2])) {
        x[, 2] <- as.character(x[, 2])
        FUN <- as.factor
    } else {
        FUN <- match.fun(paste0("as.", mode(x[, 2])))
    }	
	
    evnt <- new.env(hash = TRUE, size = nrow(x), 
        parent = emptyenv())
    outmode <- match.fun(paste0("as.", mode.out))
    apply(x, 1, function(col) {
        assign(col[1], outmode(col[2]), envir = evnt)
    })
	
    class(evnt) <- c("qdap_hash", "evir", class(evnt))
    attributes(evnt)[["mode"]] <- FUN
	evnt
} 


hash_look_e <- function(terms, envir, missing = NA) {
	
    hits <- which(!is.na(match(terms, names(as.list(envir)))))
    x <- rep(ifelse(is.null(missing), NA, missing), length(terms))
	
    x[hits] <- recoder(terms[hits], envr = envir)

    if (is.null(missing)) { 
    	keeps <- which(is.na(x))
        x[keeps] <- terms[keeps]
        x
    }   
    x
	
}


## Helper function
recoder <- function(x, envr){                               
    x <- as.character(x) #turn the numbers to character                                                        
    unlist(lapply(x, get, envir = envr))                      
}  
