# Utility functions for RcppOctave
# 
# Author: Renaud Gaujoux
# Creation: 21 Nov 2011
###############################################################################

#' @import stringr
NULL

#' Compare Lists or Environments
#' 
#' This function compares two lists or environments.
#' It is useful for comparing results obtained in R and Octave.
#'  
#' @param x a \code{list} or an \code{environment}
#' @param y a \code{list} or an \code{environment}
#' @param msg a character string used (if not missing) in a message that is 
#' printed before the comparison. It is useful for separating multiple 
#' sequential comparisons.
#' 
#' @return No value is returned, but prints out:
#' \itemize{
#' \item the element/variable names of each input list or environment, 
#' \item the result of the comparison of the elements in \code{x} and the 
#' corresponding element in \code{y} -- if present.
#' }
#' 
#' 
#' @export
#' @examples
#' 
#' X <- matrix(1:64, 8)
#' ref <- svd(X)
#' res <- .O$svd(X, argout=3)
#' 
#' check.equal(ref, res, "R and Octave function 'svd'")
#' 
check.equal <- function(x, y, msg){
	
	l1 <- x; l2 <- y; 
	if( !missing(msg) ) message('check.equal: ', msg)
	if( is.environment(l1) ) l1 <- as.list(l1)
	if( is.environment(l2) ) l2 <- as.list(l2)
	
	cat("x: ", names(l1), "\n", sep=' ')
	cat("y: ", names(l2), "\n", sep=' ')
	na <- names(l1)[names(l1) %in% names(l2)]
	sapply(names(l1), function(n){
				cat(n, ':')
				ref <- l1[[n]]
				x <- l2[[n]]
				
				if( is.matrix(ref) && ncol(ref) == nrow(ref) && ncol(ref) == 1 && is.null(dim(x)) )
					ref <- drop(ref)
				if( is.numeric(ref) ){
					storage.mode(ref) <- 'double'
					if( is.matrix(x) && any(dim(x)==1) ){
						cat('*')
						x <- as.numeric(x)
					}
				}
				if( is.numeric(x) ){
					storage.mode(x) <- 'double'
					if( is.matrix(ref) && any(dim(ref)==1) ){
						cat('*')
						ref <- as.numeric(ref)
					}
				}
				
				if( is.null(x) ){
					if( is.null(ref) || all(dim(ref) == c(0,0)) )
						cat("OK0 | ")
					else cat("NULL | ")
				}
				else if( identical(ref, x) )
					cat("OK | ")
				else{
					e <- all.equal(ref, x)
					if( is.logical(e) && e )
						cat("OK2 (", sum(abs(ref-x)),")| ")
					else
						cat("ERROR | ")
				}
				
			})
	cat("\n")
	invisible()
}

packageName <- function() 'RcppOctave' #if( pkgmaker::testRversion("> 2.15.3") ) utils::packageName else pkgmaker:::packageName

# split system PATH into a character vector
str_syspath <- function(path = Sys.getenv('PATH'), clean = FALSE){
    path <- strsplit(path, .Platform$path.sep)[[1]]
    if( clean ) path <- path[nzchar(path)]
    path
}

# prototype object to manage system PATH
Sys.path <- local({
    .init_state <- NULL # initial PATH value
    .commits <- list() # successive changes in the PATH
    
    .get <- function() Sys.getenv('PATH')
    .set <- function(x, clean = TRUE){
        if( clean ){
            sep <- .Platform$path.sep
            x <- gsub(sprintf("%s[ %s]+", sep, sep), '', x)
        }
        Sys.setenv(PATH = x)
    }
    
    list(
        get = .get
        , set = .set
        , init = function(){ # store initial state
            .init_state <<- .get()
        }
        , append = function(x){
            .set(paste(c(.get(), x), collapse = .Platform$path.sep))
        }
        , prepend = function(){
            .set(paste(c(x, .get()), collapse = .Platform$path.sep))
        }
        , rm = function(){
            
        }
        , commit = function(){
            init <- str_syspath(.init_state)
            cur <- str_syspath(.get())
            .commits <<- c(.commits, list(setdiff(cur, init)))
        }
        , revert = function(msg = NULL){
            # no message in non verbose mode
            if( !getOption('verbose') ) msg <- NULL
            
            if( !length(.commits) ) return()
            addon <- tail(.commits, 1L)[[1L]]
            if( !is.null(msg) ) message(msg, "... ", appendLF = FALSE)
            cur <- str_syspath(.get())
            p <- paste(setdiff(cur, addon), collapse = .Platform$path.sep)
            .set(p, clean = TRUE)
            if( !is.null(msg) ) message("OK")
            .commits <<- .commits[-length(.commits)]
        }
    )
            
})

# wrapper call to system (Linux) or shell (Windows) to fix an issue in
# shell when intern=TRUE and mustWork=TRUE
system_call <- function(...){
    if( .Platform$OS.type == 'windows' ){
        system <- getFunction('shell', where = 'package:base')
        res <- system(..., intern = TRUE, mustWork = TRUE)
        if( !is.null(st <- attr(res, 'status')) && st != 0 ){
            stop(paste(res, collapse = "\n  ")) 
        }
        res
    }else base::system(..., intern = TRUE)
	
}

is.Mac <- function(check.gui=FALSE){
	is.mac <- (length(grep("darwin", R.version$platform)) > 0)
	# return TRUE is running on Mac (adn optionally through GUI)
	is.mac && (!check.gui || .Platform$GUI == 'AQUA')
}

#' @importFrom utils head
file.first.path <- function(dir, ...){
    f <- file.path(...)
    sapply(f, function(x){
        x <- file.path(dir, x)
        x0 <- x[file.exists(x)]
        if( length(x0) ) head(x0, 1) else as.character(NA)
    })
    
}

