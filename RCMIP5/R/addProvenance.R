#' Add provenance information to a cmip5data object
#' 
#' It's important to track data provenance, the steps taken to produce a particular
#' dataset. Each operation in the RCMIP5 package adds provenance information via this function.
#'
#' @param x A \code{\link{cmip5data}} object
#' @param msg Either a string (message to be added to the provenance) or
#' a cmip5data object, in which case the provenance of this latter object is appended
#' to that of \code{x} (i.e., their histories are merged)
#' @param verbose logical. Print info as we go?
#' @return The original object, with an updated provenance containing:
#'  \item{timestamp}{Date and time entry was added}
#'  \item{caller}{The function that added this entry, including its parameter values}
#'  \item{message}{Description of action(s) taken}
#'  \item{dim}{Data dimensions when this entry was added}
#'  \item{digest}{Hash of the data when this entry was added; see \code{\link{digest}}}
#' @details We want to track computational steps applied to a particular
#' \code{\link{cmip5data}} object, for reproducibility and user debugging.
#' This function logs information from the caller to a 'provenance' data structure.
#' @note This is an internal RCMIP5 function and not exported.
#' @keywords internal
addProvenance <- function(x, msg, verbose=FALSE) {
    
    # Sanity checks
    stopifnot(class(x)=="cmip5data")
    stopifnot(class(msg) %in% c("character", "NULL", "cmip5data"))
    stopifnot(length(verbose)==1 & is.logical(verbose))
    
    # Get package version number, allowing that there might not be one
    pkgv <- "???"
    try({
       pkgv <- packageVersion("RCMIP5") 
    }, silent=T)
    
    if(is.null(x$provenance)) { # create a new provenance
        if(verbose) cat("Creating new provenance\n")
        x$provenance <- data.frame(
            timestamp=Sys.time(),
            caller="addProvenance",
            message=paste("RCMIP5", pkgv, "under", R.version.string),
            dim="",
            digest="",
            stringsAsFactors=F
            )
    }
    
    # Get calling function's call (its name and parameters)
    parentcall <- "<parent call unavailable>"
    try({
        parentcall <- match.call(definition=sys.function(-1), call=sys.call(-1))
        parentcall <- gsub(" ", "", paste(capture.output(parentcall), collapse=""))
        parentcall <- gsub("\\\"", "'", parentcall)
    }, silent=TRUE)
 
    # Add to the provenance data structure. Two cases: msg is a string containing
    # actual message; or it's another cmip5 data object, in which case we want to 
    # append its provenance to that of x.
    nr <- nrow(x$provenance) + 1
    x$provenance[nr, "timestamp"] <- Sys.time()
    if(class(msg) == "character") {
        # Trim multiple spaces (happens with custom FUNs in makeStat functions)
        # from http://stackoverflow.com/questions/14737266/removing-multiple-spaces-and-trailing-spaces-using-gsub
        msg <- gsub("^ *|(?<= ) | *$", "", msg, perl=T)
        # Remove artifacts for prettier code in message
        msg <- gsub("{;", "{", msg, fixed=T)
        msg <- gsub("; }", " }", msg, fixed=T)
        
        if(verbose) cat("Adding message to provenance")
        x$provenance[nr, "caller"] <- parentcall
        x$provenance[nr, "message"] <- msg
        x$provenance[nr, "dim"] <- paste(dim(x$val), collapse=",")
        dg <- "<digest unavailable>"
        try({
            dg <- digest::digest(x$val)
            } , silent=TRUE)
        x$provenance[nr, "digest"] <- dg        
    } else {
        if(verbose) cat("Appending provenances")
        x$provenance[nr, "caller"] <- "addProvenance"
        x$provenance[nr, "message"] <- "Merged (*) provenance follows:"
        yp <- msg$provenance[-1,]
        yp$message <- paste("*", yp$message)
        x$provenance <- rbind(x$provenance, yp)
    }
     
    x
} # addProvenance
