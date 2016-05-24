list_Weka_interfaces <-
function()
{
    capitalize <- function(s)
        sprintf("%s%s", toupper(substring(s, 1L, 1L)), substring(s, 2L))

    kinds <- lapply(mget(objects(Weka_interfaces), Weka_interfaces),
                    "[[", "kind")
    o <- split(names(kinds), sapply(kinds, "[", 1L))
    names(o) <-
        capitalize(sprintf("%ss",
                           sub(".*_([^_]+)_interface", "\\1", names(o))))
    writeLines(formatDL(names(o),
                        sapply(o, paste, collapse = ", "),
                        style = "list"))
    invisible(o)
}

print.R_Weka_interface <-
function(x, ...)
{
    ## We can only get the name of the interface function if print() is
    ## called explicitly.
    fname <- if(is.name(s <- substitute(x)))
        as.character(s)
    else
        ""
    
    name <- get_Java_class(x)

    if(is.function(init <- attr(x, "meta")$init))
        init()
    
    ## Seems that not all Weka learners have a globalInfo() method
    ## (e.g., Cobweb for Weka 3.4.6), so be careful.
    o <- .jnew(name)
    if(.has_method(o, "globalInfo")) {
        writeLines(c(strwrap(gettextf("An R interface to Weka class '%s', which has information",
                              as_qualified_name(name))),
                     "",
                     strwrap(.jcall(o, "S", "globalInfo"),
                             indent = 2L, exdent = 2L)))
        if(.has_method(o, "getTechnicalInformation"))
            writeLines(c("",
                         "  BibTeX:",
                         "",
                         format(get_technical_information(o),
                                offset = 2L)))
    }
    else
        writeLines(gettextf("An R interface to Weka class '%s'.",
                            as_qualified_name(name)))

    writeLines(c("",
                 gettext("Argument list:"),
                 {
                     ax <- deparse(args(x))
                     strwrap(sub("^function *", fname, ax[-length(ax)]),
                             indent = 2L, exdent = 2L)
                 }))
    
    classes <- get_R_classes_returned(x)
    if(length(classes))
        writeLines(c("",
                     gettext("Returns objects inheriting from classes:"),
                     paste(c(" ", classes), collapse = " ")))
    
    invisible(x)
}
