arrange.single.author <-
function(y)
{
    if( grepl( ",", y) ) {
        y <- sub( "^([^,]+)[[:space:]]*,[[:space:]]*(.*?)$", "\\2 \\1", y , perl = TRUE )
    }
    rx <-  "^[{](.*)[}]$"
    rx2 <- "^([^]]*)[{]([^]]*)[}]$"
    if( grepl( rx, y ) ) {
        person( sub( rx, "\\1", y ) )
    } else if( grepl( rx2, y ) ) {
        person( 
        sub( rx2, "\\1", y ), 
        sub( rx2, "\\2", y )
        )
    } else {
        as.person( y )
    }
} 

arrange.authors <-
function( x )
{
    rx <- "[[:space:]]+and[[:space:]]+"
    authors <- lapply( strsplit( x, rx )[[1]], arrange.single.author )
    as.personList( authors )
}

make.bib.entry <-
function( x )
{
    type <- attr( x, "entry" )
    key  <- attr( x, "key" )
		
    y <- as.list( x )
    names(y) <- tolower( names(y) )
		
    if( "author" %in% names(y) ){
        y[["author"]] <- arrange.authors( y[["author"]] )
    }
    if( "editor" %in% names(y) ){
        y[["editor"]] <- arrange.authors( y[["editor"]] )
    }
		
    tryCatch(  
    bibentry( bibtype = type, key = key, other = y ), 
    error = function(e){
        message( sprintf( "ignoring entry '%s' (line %d) because :\n\t%s\n", 
                         key, 
                         attr(x, "srcref")[1], 
                         conditionMessage( e ) ) )
        NULL
    } )
}

make.citation.list <-
function( x, header, footer)
{
    rval <- list()
    for( i in seq_along(x) ){
        if( !is.null(x[[i]] ) )
            rval <- c( rval, x[[i]] )
    }
    class(rval) <- c( "bibentry" )
    rval
}

findBibFile <-
function(package) {
    if( package %in% c("base", "datasets", "graphics", "grDevices", 
                       "methods", "stats", "stats4", "tools", "utils" )
       ) {
        system.file( "bib", sprintf( "%s.bib", package ), package = "bibtex" )
    } else {
        attempt <- system.file( "REFERENCES.bib", package = package )
        if( !nzchar(attempt) ){
            stop( sprintf( "no bibtex database for package '%s'", package ) ) 
        }
        attempt
    }
}


read.bib <-
function(file = findBibFile(package) , 
         package = "bibtex", 
         encoding = "unknown",
         header = if( length(preamble) ) paste( preamble, sep = "\n" ) else "", 
         footer = "" )
{
    if( !is.character( file ) ){
        stop( "'read.bib' only supports reading from files, 'file' should be a character vector of length one" )
    }
    srcfile <- switch( encoding, 
                      "unknown" = srcfile( file ), 
                      srcfile( file, encoding = encoding ) )
    out <- .External( "do_read_bib", file = file, 
                     encoding = encoding, srcfile = srcfile )
    keys <- lapply(out, function(x) attr(x, 'key'))
    at  <- attributes(out)
    if((typeof(out) != "integer") || (getRversion() < "3.0.0"))
        out <- lapply( out, make.bib.entry )
    else
        out <- list()
    preamble <- at[["preamble"]]
	
    out <- make.citation.list( out, header, footer )
    attr( out, "strings") <- at[["strings"]]
    names(out) <- keys
    out
}

#' Generate a Bibtex File from Package Citations
#'
#' Generates a Bibtex file from a list of packages or all the installed packages.
#' It is useful for adding relevant citations in Sweave documents.
#'
#' @param entry a \code{\link{bibentry}} object or a character vector of package
#' names. If \code{NULL}, then the list of all installed packages is used.
#' @param file output Bibtex file.
#' @param verbose a logical to toggle verbosity.
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
#' @export
#' @examples
#'
#' write.bib(c('bibtex', 'utils', 'tools'), file='references')
#' bibs <- read.bib('references.bib')
#' write.bib(bibs, 'references2.bib')
#' tools::md5sum(c('references.bib', 'references2.bib'))
#' md5[1] == md5[2]
#'
write.bib <-
function(entry, file="Rpackages.bib", append = FALSE, verbose = TRUE)
{
    bibs <-
    if( inherits(entry, "bibentry") )    entry
    else if( is.character(entry) ){
        if( length(entry) == 0 ){
            if( verbose ) message("Empty package list: nothing to be done.")
            return(invisible())
        }
        pkgs <- entry
        if( is.null(pkgs) ) ## use all installed packages
            pkgs <- unique(installed.packages()[,1])
        bibs <- sapply(pkgs, function(x) try(citation(x)), simplify=FALSE)
        #bibs <- lapply(pkgs, function(x) try(toBibtex(citation(x))))
        n.installed <- length(bibs)

        ## omit failed citation calls
        ok <- sapply(bibs, inherits, "bibentry")
        pkgs <- pkgs[ok]
        bibs <- bibs[ok]
        n.converted <- sum(ok)

        ## add bibtex keys to each entry
        pkgs <- lapply(seq_along(pkgs), function(i) if(length(bibs[[i]]) > 1)
                        paste(pkgs[i], 1:length(bibs[[i]]), sep = "") else pkgs[i])
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
        bibs <- mapply(function(b,k){ b$key <- k; b}, bibs, pkgs, SIMPLIFY=FALSE)
        bibs <- do.call("c", bibs)

        if(verbose) message("Converted ", n.converted, " of ", n.installed, " package citations to BibTeX")
        bibs
    } else
        stop("Invalid argument 'entry': expected a bibentry object or a character vector of package names.")

    if( length(bibs) == 0 ){
        if( verbose ) message("Empty bibentry list: nothing to be done.")
        return(invisible())
    }

    ## write everything to a single .bib file
    if( is.null(file) )
        file <- stdout()
    else if( is.character(file) ){
        if( !grepl("\\.bib$", file) ) # add .bib extension if necessary
        file <- paste(file, '.bib', sep='')
    }

    fh <- file(file, open = if(append) "a+" else "w+" )
    on.exit( if( isOpen(fh) ) close(fh) )
    if( verbose ) message("Writing ", length(bibs) , " Bibtex entries ... ", appendLF=FALSE)
    writeLines(toBibtex(bibs), fh)
    #writeLines(do.call("c", lapply(bibs, as.character)), fh)
    if(verbose) message("OK\nResults written to file '", file, "'")

    ## return Bibtex items invisibly
    invisible(bibs)
} 

