#' Wrapper for Sweave and pdf
#' 
#' This class provides the basis to build dynamic reports using LaTeX and
#' resource from a datamart. You can
#' create a report with \code{swvreport}, which takes a sweave file name. 
#' The generic method \code{put} 
#' can then be used to actually produce the report in pdf format.
#' 
#' @seealso \code{\link{swvreport}}
#' @examples
#' getSlots("SweaveReport")
#'
#' @name SweaveReport-class
#' @rdname SweaveReport-class
#' @exportClass SweaveReport
#' @author Karsten Weinert \email{k.weinert@@gmx.net}
setClass(
    Class="SweaveReport", 
    representation=representation(tpl="character"),
    contains="Target",
    validity=function(object) {
        if(tolower(strtail(object@name, 4))!=".pdf") stop("Name must end with '.pdf'")
        if(!file.exists(object@tpl)) stop("invalid tpl argument, name of existing file expected.")
    }
)

#' Constructor for SweaveReport objects
#'
#' see class SweaveReport for details.
#'
#' @param tpl         path to markdown template file
#' @param name        name of the Report, default ''
#' @param verbose     diagnostic messages T/F
#' @param clss        class of the constructed object, default 'SweaveReport'
#' @param ...         additional arguments, currently unused.
#'
#' @return generic
#' @export
#' @rdname SweaveReport-class
swvreport <- function(tpl, name=NULL, clss="SweaveReport", verbose=getOption("verbose"), ...) {
  if(is.null(name)) name <- paste(head(unlist(strsplit(basename(tpl), "\\.")),-1), ".pdf", sep="")
  # instantiate
  new(clss, name=name, tpl=tpl)
}

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,SweaveReport,DirectoryLocation-method
setMethod(
  f="put",
  signature=c(target="SweaveReport", where="DirectoryLocation"),
  definition=function(target, where, verbose=TRUE, ...) {
    if(verbose) cat("sweaving in folder '", as.character(where),"'.\n")
    old_wd <- getwd(); on.exit(setwd(old_wd))
    setwd(as.character(where))
    Sweave(target@tpl, ...)
    if(verbose) cat("working @ '", getwd(),"'.\n")
    
    setwd(as.character(where))
    tex <- paste(head(unlist(strsplit(basename(target@tpl), "\\.")),-1), ".tex", sep="")
    if(verbose) cat("texifying '", tex,"'.\n")
    if(!file.exists(tex)) stop("Could not find Sweave output ", tex, ". Did Sweave fail?")
    cmd <- paste("texify -p -b -q", basename(tex))
    system(cmd)
    res_from <- paste(strhead(tex, -3), "pdf", sep="")
    if(!file.exists(res_from)) stop("Could not find PDF output ", res_from, ". Did texify fail?")
    res_to <- file.path(as.character(where), target@name)
    file.rename(res_from, res_to)
    return(res_to)
  }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,SweaveReport,missing-method
setMethod(
  f="put",
  signature=c(target="SweaveReport", where="missing"),
  definition=function(target, where, verbose=TRUE, ...) {
    if(verbose) cat("assigning folder '", dirname(target@tpl),"'.\n")
    put(target=target, where=dirloc(path=dirname(target@tpl)), ...)
  }
)
