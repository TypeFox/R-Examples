WOW <-
function(x)
{
    ## Weka Option Wizard.

    name <- get_Java_class(x)
    if(is.null(name))
        stop("Argument 'x' does not specify a known Weka method.")

    if(inherits(x, "R_Weka_interface")) {
        if(is.function(init <- attr(x, "meta")$init))
            init()
    }
    
    x <- .jnew(name)

    names <- descriptions <- synopses <- character()
    lengths <- integer()

    if(.has_method(x, "listOptions")) {
        opt <- .jcall(x, "Ljava/util/Enumeration;", "listOptions")
        while(.jcall(opt, "Z", "hasMoreElements")) {
            o <- .jcall(opt, "Ljava/lang/Object;", "nextElement")
            ## In fact, o is now a weka.core.Option object.
            names <- c(names, .jcall(o, "S", "name"))
            descriptions <-
                c(descriptions, .jcall(o, "S", "description"))
            synopses <- c(synopses, .jcall(o, "S", "synopsis"))
            lengths <- c(lengths, .jcall(o, "I", "numArguments"))
        }
    }

    `class<-`(list(Name = names,
                   Length = lengths,
                   Description = descriptions,
                   Synopsis = synopses),
              "WOW")
    
}
    
print.WOW <-
function(x, ...)
{
    if(length(x$Name)) {
        ## Note that we get nothing in case the Weka class has no
        ## listOptions() method.
        out <- mapply(formatDL,
                      x$Synopsis,
                      gsub("\t", " ", x$Description),
                      indent = 8L)
        if(any(ind <- (x$Length > 0L)))
            out[ind] <- gettextf("%s\n\tNumber of arguments: %d.",
                                 out[ind], x$Length[ind])
        writeLines(out)
    }
    else
        writeLines(gettext("No available information on options."))
    invisible(x)
}
