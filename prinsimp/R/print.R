### Print methods for the simpart class

print.simpart <- function(x, ...) {
    d <- nrow(x$simple)
    modeldim <- ncol(x$model)
    simpledim <- ncol(x$simple)
    
    if (!is.null(x$call))
        cat(paste0('\nCall:\n', deparse(x$call), '\n\n'))

    cat(paste0("Simplicity measure: ",
               if (!is.null(x$measure)) {
                   if (x$measure %in% c('first', 'second')) {
                       paste(x$measure, 'divided differences')
                   } else x$measure
               } else 'unknown',
               "\n\n"))
    
    cat(paste0('Partition simplicity (' , simpledim, ' simple basis):\n'))
    simplicity <- c(x$simplicity$model, x$simplicity$simple)

    names(simplicity) <- paste(rep(c('model', 'simple'), c(modeldim, simpledim)),
                               c(seq_len(modeldim), seq_len(simpledim)))
    print(simplicity, ...)

    cat(paste0('\nFull space simplicity:\n'))
    full_simplicity <- x$simplicity$full
    names(full_simplicity) <- paste('full', seq_len(d))
    print(full_simplicity, ...)

    invisible(x)
}


summary.simpart <- function (object, loadings = FALSE, ...) {
    object$print.loadings <- loadings
    class(object) <- "summary.simpart"
    object
}


print.summary.simpart <- function (x, digits = 3, loadings = x$print.loadings, ...) {
    d <- nrow(x$simple)
    modeldim <- ncol(x$model)
    simpledim <- ncol(x$simple)

    cat(paste0('Simple partition (',
               if (!is.null(x$measure)) {
                   if (x$measure %in% c('first', 'second')) {
                       paste(x$measure, 'divided differences')
                   } else x$measure
               } else 'unknown',
               '): ', simpledim, ' simple basis\n\n'))
    
    nm <- matrix(c(rep(c(' model', 'simple'),
                       c(d-simpledim, simpledim)),
                   c(seq_len(d-simpledim), seq_len(simpledim))),
                 ncol=2)
    names <- apply(nm, 1, function(r) paste(r[1], r[2]))

    width <- max(nchar(names, type='w'),
                 nchar(format(unlist(x$simplicity[1:2]), digits=digits), type='w'))

    
    vals <- rbind(`Simplicity` = format(unlist(x$simplicity[1:2]),
                                        digits = digits, width = width),
                  `%-var explained` = format_varperc(unlist(x$varperc),
                                                     digits = digits, width = width),
                  `Cumulative %-var` = format_varperc(c(cumsum(x$varperc$model),
                                                        cumsum(x$varperc$simple)),
                                               digits = digits, width = width))
    colnames(vals) <- paste(rep(c('model', 'simple'), c(modeldim, simpledim)),
                            c(seq_len(modeldim), seq_len(simpledim)))
    print(vals, digits = digits, quote = FALSE, ...)

    if (loadings) {
        basis <- as.data.frame(cbind(x$model, x$simple))
        colnames(basis) <- names
        cat("\nLoadings:\n")
        cx <- format(round(basis, digits = digits))
        print(cx, quote = FALSE, ...)
    }

    invisible(x)
}


### Formats the percentage variance explained (i.e., in range 0-100)
### to a limited width, replacing small numbers with the '<.1' label.
###
### Returns the character vector of formatted labels
format_varperc <- function(x, digits = 3, width = NULL, ...) {
    xf <- prettyNum(x, digits = digits, width = width, ...)
    xf[x < 0.1] <- paste0(paste0(rep(' ', max(0, width-3)),
                                                     collapse=''),
                                              '<.1', collapse='')
    xf
}
