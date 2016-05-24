#' Buildable markdown report
#' 
#' This class provides the basis for building dynamic reports using
#' markdown text and resource from a \link{datamart}. You can
#' create a report with \code{mdreport}, which takes a template file name,
#' and a list of variables as arguments. The generic method \code{put} 
#' can then be used to actually produce the report at various locations
#' (directory, memory, blogging site).
#'
#' \code{strsubst} is a simple
#' templating mechanism inspired from Python (PEP-0292).
#' Variables in the template are marked by a preceeding
#' dollar sign and get replaced with the value
#' of the corresponding variables passed to \code{strsubst}.
#'
#' @seealso \code{\link{mdreport}}, \code{\link{swvreport}}
#'
#' @examples
#' getSlots("MdReport")
#'
#' @name MdReport-class
#' @rdname MdReport-class
#' @exportClass MdReport
setClass(
    Class="MdReport", 
    representation=representation(tpl="character", vars="list", xdata="Xdata"),
    contains="Target"
)

#' Constructor for MdReport objects
#'
#' see class MdReport for details.
#'
#' @param tpl         path to markdown template file
#' @param name        name of the Report, default ''
#' @param xdata       instance of Xdata class
#' @param verbose     diagnostic messages T/F
#' @param clss        class of the object, default 'MdReport'
#' @param ...         number of targets
#'
#' @return generic
#' @export
#' @rdname MdReport-class
mdreport <- function(tpl, name="", xdata=NULL, clss="MdReport", verbose=getOption("verbose"), ...) {
    if(is.null(xdata)) xdata <- new("EmptySet")
    # handle subtargets
    vars <- list(...)
    idx <- which(names(vars)=="")
    for (i in idx) {
        v <- vars[[i]]
        if(is.character(v@name) && v@name!="") names(vars)[[i]] <- v@name
    }
    if(verbose) {
    if(length(vars)==0) 
        cat("no subtargets.\n")
    else
        cat("subtargets: '", paste(names(vars), collapse="', '"), "'.\n")
    }

    # tpl filename?
    if(length(tpl)==1 && file.exists(tpl)) tpl <- readLines(tpl)
    if(!is.character(tpl)) stop("invalid tpl argument, name of existing file or character vector expected.")

    # interpret the first lines as meta data and convert them to subtargets  
    vname <- ""
    body_start <- 1
    for (l in tpl) {
        body_start <- body_start + 1
        if(grepl("^[\\s]*$", l, perl=TRUE)) break
        m <- strparse("^[\\s]+(?<vval>.*)", l)[[1]] # continuation line?
        if(!is.null(m)) {
            if(vname=="") stop("invalid header.")
            vars[[vname]] <- paste(vars[[vname]], m[["vval"]])
            next
        }
        m <- strparse("^(?<vname>[^:\\s]+)[\\s]*:[\\s]*(?<vval>.*)$", l)[[1]] # meta?
        if(!is.null(m)) {
            vname <- strcap(m[["vname"]])
            vars[[vname]] <- m[["vval"]]
        } else if(body_start==2) { # no meta
            if(verbose) cat("no meta data found.\n")
            body_start <- 1
            break
        }
    }
    if(verbose) {
    cat(sprintf("processed %d header lines.\n", body_start-1))
    if(length(vars)==0) 
        cat("still no subtargets.\n")
    else
        cat("subtargets after meta: '", paste(names(vars), collapse="', '"), "'.\n")
    }

    # body
    if(body_start>1) tpl <- tail(tpl, -body_start+1)
    tpl <- paste(tpl, collapse="\n")

    # instantiate
    new(clss, name=name, tpl=tpl, vars=vars, xdata=xdata)
}

#' @param draft    Parameter to MdReport, is draft? Logical.
#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,MdReport,Location-method
setMethod(
    f="put",
    signature=c(target="MdReport", where="Location"),
    definition=function(target, where, draft=TRUE, overwrite=TRUE, ...) {
        pat <- "\\$\\((?<name>[^\\)]+)\\)"
        mem <- new("MemoryLocation")
        tgts <- unique(unlist(strparse(pat, target@tpl)))
        vars <- target@vars
        for (tg in tgts) {
            # browser()
            tgv <- vars[[tg]]
            if(inherits(tgv, "Target")) vars[[tg]] <- put(tgv, mem)
            # if(is.null(tgv)) 
            # vars[[tg]] <- as.character(put(target=v, where=mem))
        }

        mdlines <- strsubst(target@tpl, vars)
        html <- markdown::markdownToHTML(
            text=mdlines,
            options=c("safelink", "escape", "use_xhtml", "fragment_only", "smartypants"),
            extensions=NULL
        )
        
        subj <- vars[["Subject"]]
        if(inherits(subj, "Target")) subj <- put(target=subj, where=mem)
        if(!is.character(subj)) subj <- "Untitled"
        
        keyw <- target@vars[["Keywords"]]
        if(inherits(keyw, "Target")) keyw <- put(target=keyw, where=mem)
        if(!is.character(keyw)) keyw <- ""
        
        put(blogpost(name=target@name, subject=subj, content=html, label=keyw, draft=draft, overwrite=overwrite, ...), where)
    }
) 


