unf <- function(x, version = 6, ...){
    if(is.matrix(x))
        x <- as.data.frame(x)
    if((is.data.frame(x) | is.list(x)) & length(x)==1)
        x <- x[[1]]
    if(is.data.frame(x) | is.list(x)){
        locale <- Sys.getlocale(category="LC_COLLATE")
        Sys.setlocale(category="LC_COLLATE", "C")
        if(version==3){
            vars <- sapply(x, function(i) unf3(i, ...)$unf)
            out <- unf3(sort(vars), ...)
        } else if(version==4){
            vars <- sapply(x, function(i) unf4(i, ...)$unf)
            out <- unf4(sort(vars), ...)
        } else if(version==4.1){
            vars <- sapply(x, function(i) unf4(i, version = 4.1, ...)$unf)
            out <- unf4(sort(vars), version = 4.1, ...)
        } else if(version==5){
            vars <- sapply(x, function(i) unf5(i, ...)$unf)
            out <- unf5(sort(vars), ...)
        } else if(version==6){
            vars <- sapply(x, function(i) unf6(i, ...)$unf)
            out <- unf6(sort(vars), ...)
        } else {
            stop("Unrecognized UNF version: must be 3, 4, 4.1, 5, or 6.")
        }
        Sys.setlocale(category="LC_COLLATE", locale)
        out$variables <- vars
        return(out)
    } else {
        if(version==3){
            out <- unf3(x, ...)
        } else if(version==4){
            out <- unf4(x, ...)
        } else if(version==4.1){
            out <- unf4(x, version=4.1, ...)
        } else if(version==5){
            out <- unf5(x, ...)
        } else if(version==6){
            out <- unf6(x, ...)
        } else {
            stop("Unrecognized UNF version: must be 3, 4, 4.1, 5, or 6.")
        }
        return(out)
    }
}

print.UNF <- function(x, ...){
    if('formatted' %in% names(x)) {
        out <- x$formatted
    } else {
        if(is.null(attr(x,'version'))) {
            out <- paste0('UNF:', x$unf)
        } else {
            if(attr(x, 'version')==6) {
                digits <- ifelse(!is.null(attr(x, "digits")), attr(x, "digits"), 7)
                characters <- ifelse(!is.null(attr(x, "characters")), attr(x, "characters"), 128)
                truncation <- ifelse(!is.null(attr(x, "truncation")), attr(x, "truncation"), 128)
                header <- paste(if(digits != 7) paste0("N", digits) else NULL,
                                if(characters != 128) paste0("X", characters) else NULL,
                                if(truncation != 128) paste0("H", truncation) else NULL,
                                sep = ",", collapse="")
                header <- ifelse(length(header), gsub("^[[:punct:]]+", "", header), "")
                header <- ifelse(length(header), gsub("[[:punct:]]+$", "", header), "")
                header <- ifelse(length(header), gsub("[[:punct:]]{2}", ",", header), "")
                out <- paste0('UNF6:', ifelse(header == "", x$unf, paste0(header,':', x$unf)))
            } else if(attr(x, 'version') %in% c(3,4,4.1,5)) {
                out <- paste0('UNF',version,':',
                    if((!is.null(attr(x,'digits')) & attr(x,'digits')!=7) |
                        (!is.null(attr(x,'characters')) & attr(x,'characters')!=128)) {
                        paste0(paste(attr(x,'digits'), attr(x,'characters'), sep=','), ':', x$unf)
                    } else {
                        x$unf
                    })
            } else {
                out <- paste0('UNF', attr(x, 'version'), ':', x$unf)
            }
        }
    }
    cat(out, '\n')
    invisible(x)
}
