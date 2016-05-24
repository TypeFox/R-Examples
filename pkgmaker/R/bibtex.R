# Bibtex related functions
# 
# Author: Renaud Gaujoux
# Creation: 25 Apr 2012
###############################################################################


#' Generate a Bibtex File from Package Citations
#' 
#' Generates a Bibtex file from a list of packages or all the installed packages.
#' It is useful for adding relevant citations in Sweave documents.
#' 
#' Multiple citations are handled by adding a numeric suffix to the Bibtex key 
#' (other than the first/main citation) as \code{"<pkgname>\%i"} (e.g. pkg, pkg2, pkg3).
#' 
#' This function has now been integrated by Romain François in the bibtex package.
#'
#' @encoding utf8
#'
#' @param entry a \code{\link{bibentry}} object or a character vector of package 
#' names. If \code{NULL}, then the list of all installed packages is used.
#' @param file output Bibtex file. It can be specified as a filename (as a single 
#' character string), NULL for \code{stdout}, or a \code{link{connection}} object. 
#' If \code{file} is a character string, an extension '.bib' is appended if not 
#' already present.
#' @param prefix character string to prepend to the generated packages' Bibtex key.
#' @param append a logical that indicates that the Bibtex entries should be added
#' to the file. If \code{FALSE} (default), the file is overwritten.  
#' @param verbose a logical to toggle verbosity. If \code{file=NULL}, verbosity 
#' is forced off. 
#'
#' @return the list of Bibtex objects -- invisibly.
#' @author
#' Renaud Gaujoux, based on the function \code{Rpackages.bib} 
#' from Achim Zeileis (see \emph{References}).
#' 
#' @references 
#' \emph{[R] Creating bibtex file of all installed packages?}
#' Achim Zeileis. R-help mailing list. 
#' \url{https://stat.ethz.ch/pipermail/r-help/2009-December/222201.html}
#' 
#' @seealso \code{link{connection}}, \code{link{bibentry}}
#'  
#' @export
#' @examples
#' 
#' write.pkgbib(c('bibtex', 'utils', 'tools'), file='references')
#' bibs <- bibtex::read.bib('references.bib')
#' write.pkgbib(bibs, 'references2.bib')
#' md5 <- tools::md5sum(c('references.bib', 'references2.bib'))
#' md5[1] == md5[2]
#' \dontshow{ stopifnot(md5[1] == md5[2]) }
#' 
#' # write to stdout()
#' write.pkgbib(c('bibtex', 'utils', 'tools'), file=NULL)
#' 
#' # clean up 
#' unlink(c('references.bib', 'references2.bib'))
#' 
write.pkgbib <- function(entry=NULL, file="Rpackages.bib", prefix='', append = FALSE, verbose = TRUE)
{
	# special handling of file=NULL: use stdout()
	if( is.null(file) ){
		file <- stdout()
		verbose <- FALSE
	}	
	## use all installed packages if nothing is specified
	if( is.null(entry) ){ 
		if( verbose ) message("Generating Bibtex entries for all installed packages ", appendLF=FALSE)
		entry <- unique(installed.packages()[,1])
		if( verbose ) message("[", length(entry), "]")
	}
	
	bibs <- 
			if( is(entry, 'bibentry') )	entry
			else if( is.character(entry) ){
				if( length(entry) == 0 ){
					if( verbose ) message("Empty package list: nothing to be done.")
					return(invisible())
				}
				
				pkgs <- entry
				bibs <- sapply(pkgs, function(x) try(citation(x)), simplify=FALSE)
				#bibs <- lapply(pkgs, function(x) try(toBibtex(citation(x))))
				n.installed <- length(bibs)
				
				## omit failed citation calls
				ok <- sapply(bibs, is, 'bibentry')
				pkgs <- pkgs[ok]
				bibs <- bibs[ok]
				n.converted <- sum(ok)
				
				## add bibtex keys to each entry
				pkgs <- lapply(seq_along(pkgs), function(i){
							if(length(bibs[[i]]) > 1)
								paste(prefix, pkgs[i], c('', 2:length(bibs[[i]])), sep = "") 
							else paste(prefix, pkgs[i], sep='')
				})
				pkgs <- do.call("c", pkgs)
				bibs <- do.call("c", bibs)		
				# formatting function for bibtex keys:
				# names with special characters must be enclosed in {}, others not.
				as.bibkey <- function(x){
					i <- grep("[.]", x)
					if( length(i) > 0 )
						x[i] <- paste("{", x[i], "}", sep='')
					x
				}		
				#bibs <- mapply(function(b,k){ if( is.null(b$key) ) b$key <- as.bibkey(k); b}, bibs, pkgs, SIMPLIFY=FALSE)
				bibs <- mapply(function(b,k){ if( is.null(b$key) ) b$key <- k; b}, bibs, pkgs, SIMPLIFY=FALSE)
				bibs <- do.call("c", bibs)
				
				if(verbose) message("Converted ", n.converted, " of ", n.installed, " package citations to BibTeX")					
				bibs
			} else
				stop("Invalid argument `entry`: expected a bibentry object or a character vector of package names.")
	
	if( length(bibs) == 0 ){
		if( verbose ) message("Empty bibentry list: nothing to be done.")
		return(invisible())
	}
	
	## write everything to the .bib file
	not_anonymous <- !identical(file,'')
	fh <- if( is.character(file) ){
				if( not_anonymous && !grepl("\\.bib$", file) ) # add .bib extension if necessary 
					file <- paste(file, '.bib', sep='')
				fh <- file(file, open = if(append && not_anonymous) "a+" else "w+" )
				if( not_anonymous )
					on.exit( if( isOpen(fh) ) close(fh) )
				fh
			} else if( is(file, 'connection') )
				file
			else
				stop("Invalid argument `file`: expected a filename, NULL, or a connection [", class(file), "]")
	
	if( !is(fh, 'connection') )
		stop("Invalid connection: ", fh)		
	file.desc <- summary(fh)['description']
	
	if( verbose ) message(if( append ) "Adding " else "Writing ", length(bibs) , " Bibtex entries ... ", appendLF=FALSE)
	writeLines(toBibtex(bibs), fh)
	if(verbose) message("OK\nResults written to file '", file.desc, "'")
	
	## return Bibtex items invisibly
	if( !not_anonymous ) attr(bibs, 'connection') <- fh 
	invisible(bibs)
}

#' @rdname pkgmaker-defunct
write.bib <- function(...){
	.Defunct('write.pkgbib', package = 'pkgmaker')
}

#' Package References
#' 
#' Create a citation string from package specific BibTex entries, suitable to 
#' be used in Rd files.
#' The entries are looked in a file named REFERNCES.bib in the package's root 
#' directory (i.e. inst/ in development mode).
#'  
#' @param key character vector of BibTex keys
#' @param short logical that indicates if the reference should be shorten as 
#' First Author et al. if it has more than one author.
#' @return a character string containing the text formated BibTex entries
#'  
#' @export
packageReference <- function(key, short=FALSE){
	bibs <- bibtex::read.bib(file=packageReferenceFile())
	k <- sapply(bibs, function(x) x$key)
    mk <- match(key, k)
	sel <- mk[!is.na(mk)]
	if( !length(sel) ) return("")
	if( !short ){
		paste(format(bibs[sel]), collapse="\n\n")
	}else{
		sapply(bibs[sel], function(x){
					if( length(x$author$family) <= 1L ) 
						paste(x$author$family, '(', x$year, ')', sep='')				
					else{
						paste(x$author$family[[1]], ' et al. (', x$year, ')', sep='')
					}
				})
	} 
}

#' Citing Package References
#' 
#' Create a citation command from package specific BibTex entries, suitable to 
#' be used in Rd files or Latex documents.
#' The entries are looked in a file named REFERNCES.bib in the package's root 
#' directory (i.e. inst/ in development mode).
#'  
#' @param key character vector of BibTex keys
#' @param ... extra arguments passed to \code{format.bibentry}.
#' @param REFERENCES package or bibentry specification
#' @return a character string containing the text formated BibTex entries
#'  
#' @export
#' 
citecmd <- local({
	
	.init <- list(REFERENCES=NULL, KEYS=NULL) 
	.cache <- .init
	function(key, ..., REFERENCES=NULL){
		 
		# detect package name if necessary
		if( is.null(REFERENCES) ){
			# reset if explicitly passed NULL
			if( hasArg(REFERENCES) ) .cache <<- .init
			
			if( is.null(.cache$REFERENCES) ){
				pkg <- Sys.getenv('R_PACKAGE_NAME')
				if( !nchar(pkg) )
					pkg <- Sys.getenv('R_INSTALL_PKG')
				if( !nchar(pkg) )
					pkg <- Sys.getenv('MAKE_R_PACKAGE')
				if( !nchar(pkg) )
					stop("Could not identify package")
				# load REFERENCES from detected package
				.cache$REFERENCES <<- bibtex::read.bib(package=pkg)
			}
			REFERENCES <- .cache$REFERENCES
		}
		
		# load relevant Bibtex file
		REFERENCES <- if( is(REFERENCES, 'bibentry') ) REFERENCES
				else if( is.character(REFERENCES) ){
					p <- str_match(REFERENCES, "^package:(.*)")[,2]
					if( is.na(p) ) bibtex::read.bib(file=REFERENCES)
					else bibtex::read.bib(package=p)
				}else
					stop("Invalid argument `REFERENCES`: expected bibentry object or character string [", class(REFERENCES), "]")
		
		# update the cache if no keys are provided
		if( missing(key) ){
			.cache$REFERENCES <<- REFERENCES
			if( hasArg(REFERENCES) ) return(invisible(.cache$KEYS))
			else return(.cache$KEYS)
		}
		
		# check key type
		if( !is.character(key) )
			stop("Invalid argument `key`: must be a character vector.")
		
		# extract the Bibtex keys
		refkey <- sapply(REFERENCES, function(x) x$key)
		pkgs <- str_match(key, "^package:(.*)")[,2]
		nokey <- !key %in% refkey
		i_pkgs <- which(nokey && !is.na(pkgs))
		if( length(i_pkgs) > 0L ){
			# only include \cite{key} if running Sweave
			.cache$KEYS <<- unique(c(.cache$KEYS, key[i_pkgs]))
			key[i_pkgs] <- pkgs[i_pkgs] 
		}
		paste("\\cite{", key, "}", sep='')
#		if( inSweave() ) paste("\\cite{", k, "}", sep='')
#		else paste(format(REFERENCES[k %in% key], ...), collapse="\n\n")
	}
})

citecmd_pkg <- function(key, ...){
	citecmd(str_c('package:', key), ...)
}

#' Bibtex Utilities
#' 
#' \code{packageReferenceFile} returns the path to a package REFERENCES.bib file.
#' 
#' @param PACKAGE package name
#' 
#' @rdname bibtex
packageReferenceFile <- function(PACKAGE=NULL) packagePath('REFERENCES.bib', package=PACKAGE)

