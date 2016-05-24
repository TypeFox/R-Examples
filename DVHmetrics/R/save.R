saveMetric <-
function(x, file="", ...) {
    UseMethod("saveMetric")
}

saveMetric.data.frame <-
function(x, file="", ...) {
    write.table(x, file=file, row.names=FALSE, ...)
}

saveMetric.list <-
function(x, file="", ...) {
    addArgs <- list(...)

    ## determine sensible file names
    pat1st <- "([[:alnum:]]+[[:alnum:][:blank:][:punct:]]*)\\.[[:alnum:]]+$"
    pat2nd <-  "[[:alnum:]]+[[:alnum:][:blank:][:punct:]]*(\\.[[:alnum:]]+$)"

    first <- if(grepl(pat1st, file)) {
        first <- paste0(sub(pat1st, "\\1", file), "_", collapse="")
    } else {
        "metric_"
    }

    last <- if(grepl(pat2nd, file)) {
        sub(pat2nd, "\\1", file)
    } else {
        ".txt"
    }

    comps  <- gsub("\\.", "_", names(x))
    fNames <- paste0(first, comps, last)

    if(length(addArgs) > 0) {
        Map(saveMetric.data.frame, x, file=fNames, ...)
    } else {
        Map(saveMetric.data.frame, x, file=fNames)
    }

    return(invisible(NULL))
}

saveDVH <-
function(x, file="", ...) {
    addArgs <- list(...)

    ## determine sensible file names
    pat1st <- "([[:alnum:]]+[[:alnum:][:blank:][:punct:]]*)\\.[[:alnum:]]+$"
    pat2nd <-  "[[:alnum:]]+[[:alnum:][:blank:][:punct:]]*(\\.[[:alnum:]]+$)"

    first <- if(grepl(pat1st, file)) {
        first <- paste0(sub(pat1st, "\\1", file), "_", collapse="")
    } else {
        "dvh_"
    }

    last <- if(grepl(pat2nd, file)) {
        sub(pat2nd, "\\1", file)
    } else {
        ".pdf"
    }

    comps <- gsub("\\.", "_", names(x))
    file  <- paste0(first, comps, last)

    ### ... my be NULL -> cannot use Map()
    if(length(addArgs) > 0) {
        Map(ggsave, x, filename=file, ...)
    } else {
        Map(ggsave, x, filename=file)
    }

    return(invisible(NULL))
}

saveConstraint <-
function(x, ...) {
    write.table(x, row.names=FALSE, ...)
}
