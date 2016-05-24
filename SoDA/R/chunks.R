## R functions to interface to Perl routines
## The Perl code routines must have been defined
## See the comments in tests/chunk*.R

chunksAdd <- function(table, data, convert)
    stop("This function has been decommissioned to avoid hassles from CRAN caused by (protected) calls to the RSPerl package.  If you want to experiment with it, uncomment the code at the bottom of this file and in the .onLoad function.")

chunksDrop <- function(table, data, convert)
    stop("This function has been decommissioned to avoid hassles from CRAN caused by (protected) calls to the RSPerl package.  If you want to experiment with it, uncomment the code at the bottom of this file and in the .onLoad function.")


## chunksadd <- function( table = RSPerl::.PerlExpr("\\%{0};", .convert = FALSE),
##                       data = character(),
##                       convert = length(data) == 0) {
##     if(!inherits(table, "PerlHashReference"))
##       stop(gettextf(
##        "Argument table must be reference to a Perl hash object; got an object of class \"%s\"",
##                     class(table)), domain = NA)
##     args <- c(list(table), as.list(data))
##  RSPerl::.Perl("chunks_add",  .args = args, convert = convert)
## }

## chunksDrop <- function(table,
##                        data = character(),
##                        convert = length(data) == 0) {
##     if(missing(table))
##       stop("Must start with a non-empty table")
##     else if(!inherits(table, "PerlHashReference"))
##       stop(gettextf(
##        "Argument table must be reference to a Perl hash object; got an object of class \"%s\"",
##                     class(table)), domain = NA)
##     args <- c(list(table), as.list(data))
##     RSPerl::.Perl("chunks_drop",  .args = args, convert = convert)
## }
