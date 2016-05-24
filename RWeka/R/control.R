Weka_control <-
function(...)
{
    rval <- list(...)
    if((length(rval) > 0L)
       && (is.null(names(rval)) ||
           !all(nzchar(names(rval)) | sapply(rval, identical, "--"))))
        stop("All arguments must be named.")
    `class<-`(rval, "Weka_control")
}

print.Weka_control <-
function(x, ...)
{
    if(length(x) < 1L) {
        writeLines(gettext("An empty list of Weka control options."))
    } else {
        writeLines(gettextf("A list of Weka control options (%s):\n",
                            paste(as.character(x), collapse = " ")))
        ## (Note that this is quite perfect in case of control lists
        ## with recursive elements.)
        print.default(unclass(x))
    }
    invisible(x)
}

as.character.Weka_control <-
function(x, ...)
{
    if(!length(x)) return(character())

    ## <NOTE>
    ## R 2.6.0 has base::Map(), but (currently?) this has
    ## USE.NAMES = TRUE ...
    map <- function(f, ...)
        mapply(f, ..., SIMPLIFY = FALSE, USE.NAMES = FALSE)
    ## </NOTE>

    arg2char <- function(tag, val) {
        ## Works for a single (possibly empty) tag and a single value,
        ## which can be recursive, so that we can handle things like
        ##   Weka_control(K = list("RBFKernel", G = 2))
        ## <FIXME>
        tag <- if(nzchar(tag)) paste("-", tag, sep = "") else NULL
        out <- if(is.list(val)) {
            nms <- names(val)
            if(is.null(nms)) nms <- character(length(val))
            c(tag,
              paste(unlist(map(arg2char, nms, val)), collapse = " "))
        }
        else if(is.logical(val)) {
            if(val) tag
        }
        else if(inherits(val, "R_Weka_interface"))
            c(tag, get_Java_class(val))
        else
            c(tag, as.character(val))
        out
    }

    unlist(map(arg2char, names(x), x))
}

## Handlers.

make_Weka_control_handler <-
function(options, fun, ...)
{
    ## Return a function which applies 'fun' to the values of the
    ## control arguments with names in options.
    ##
    ## This is useful as e.g. for meta learners, '-W' (for giving the
    ## base learner) requires the full Java class name, but R/Weka users
    ## do not necessarily (have to) know this.  Similarly, for savers,
    ## '-c' (for giving the class index) uses Java-style counting for
    ## its attributes, starting at 0.

    ## We handle both new-style 'control' arguments to R/Weka interface
    ## functions given via RWeka_control() and old-style specifications
    ## as character vectors.  For the former, the option values can be
    ## arbitrary R objects.

    ## This is not intended to validate the given control options.

    ## Force promises:
    options; fun
    
    function(x) {
        if(inherits(x, "Weka_control")) {
            ## New-style Weka control lists:
            ##   Weka_control(W = value)
            ind <- which(names(x) %in% substring(options, 2L))
            if(any(ind)) x[ind] <- lapply(x[ind], fun, ...)
            ## Do *not* coerce to character so that we can compose
            ## several control handlers.
        }
        else {
            ## Old-style character stuff:
            ##   c("-W", "J48")
            x <- as.character(x)        # Just making sure ...
            ind <- which(x %in% options)
            if(any(ind)) x[ind + 1L] <- sapply(x[ind + 1L], fun, ...)
        }
        x
    }
}

.control_handlers <-
function(olist, flist)
{
    ## A little utility to enhance code legibility.
    if(!is.list(olist)) olist <- list(olist)
    if(!is.list(flist)) flist <- list(flist)
    if((len <- length(olist)) != length(flist))
        stop("Argument lengths differ.")
    list(control = Map(make_Weka_control_handler, olist, flist))
}
