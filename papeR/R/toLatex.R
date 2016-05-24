## functions for coercion to LaTeX "objects"

## overwrite standard generic
toLatex <- function(object, ...)
    UseMethod("toLatex")

## per default fall back to standard generic
toLatex.default <- function(object, ...)
    utils::toLatex(object, ...) # nocov

## modified version based on sanitize subroutine defined in
## R package xtable (Version  1.7-1) inside function print.xtable
##
## URL of original package: http://CRAN.R-project.org/package=xtable
## Authors of R package xtable (inlcuding print.xtable):
##   David Dahl <email: dahl@stat.tamu.edu> with contributions and
##   suggestions from many others (see source code).
##
## Licence of R package xtable: GPL-2 | GPL-3
toLatex.character <- function(object, ...) {
    result <- object
    result <- gsub("$\\sum$", "SUM", result, fixed = TRUE)
    result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
    result <- gsub("$", "\\$", result, fixed = TRUE)
    result <- gsub(">=", "$\\geq$", result, fixed = TRUE)
    result <- gsub("<=", "$\\leq$", result, fixed = TRUE)
    result <- gsub(">", "$>$", result, fixed = TRUE)
    result <- gsub("<", "$<$", result, fixed = TRUE)
    result <- gsub("|", "$|$", result, fixed = TRUE)
    result <- gsub("{", "\\{", result, fixed = TRUE)
    result <- gsub("}", "\\}", result, fixed = TRUE)
    result <- gsub("%", "\\%", result, fixed = TRUE)
    result <- gsub("&", "\\&", result, fixed = TRUE)
    result <- gsub("_", "\\_", result, fixed = TRUE)
    result <- gsub("#", "\\#", result, fixed = TRUE)
    ## a^NOT_A_NUMBER
    result <- gsub("\\^([^[:digit:]])", "\\\\verb|^|\\1", result)
    result <- gsub("\\^([[:digit:]]+)", "$^{\\1}$", result)
    result <- gsub("~", "\\~{}", result, fixed = TRUE)
    ## grep for ^2 and ^3
    result <- gsub("\u00B2", "$^2$", result, fixed = TRUE)
    result <- gsub("\u00B3", "$^3$", result, fixed = TRUE)
    result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$",
                   result, fixed = TRUE)
    result <- gsub("SUM", "$\\sum$", result, fixed = TRUE)
    return(result)
}


## modified version based on toLatex.sessionInfo from package utils
##
## Copyright (C) 1995-2013 The R Core Team
## URL: http://cran.at.r-project.org/src/base/R-3/R-3.0.1.tar.gz
## Inside archive path: /src/library/utils/R/sessionInfo.R
## Licence of R package utils: >= GPL-2
##
## with major changes and modifications by Benjamin Hofner
toLatex.sessionInfo <- function(object, pkgs = NULL, locale = FALSE,
                                base.pkgs = FALSE, other.pkgs = TRUE,
                                namespace.pkgs = FALSE, citations = TRUE,
                                citecommand = "\\citep", file = NULL,
                                append = FALSE, ...) {
    if (!is.null(pkgs)) {
        object <- sessionInfo(package = pkgs)
        if (!other.pkgs)
            warning(sQuote("other.pkgs"), " should be TRUE if ",
                    sQuote("pkgs"), " is specified.")
    }

    opkgver <- sapply(object$otherPkgs, function(x) x$Version)
    nspkgver <- sapply(object$loadedOnly, function(x) x$Version)
    key <- NULL

    if (citations) {
        bibs <- write.bib("base", file = file, append = append, verbose = FALSE)
        all_bibs <- bibs
        key <- bibs$key
    }

    z <- c("\\begin{itemize}\\raggedright",
           paste0("  \\item ", object$R.version$version.string,
                  if (citations)
                      paste0(citecommand, "{", key, "}")))

    if (locale) {
        z <- c(z, paste0("  \\item Locale: \\verb|",
                         gsub(";", "|, \\\\verb|", object$locale), "|"))
    }

    if (base.pkgs) {
        z <- c(z, strwrap(paste("\\item Base packages: ",
                                paste(sort(object$basePkgs), collapse = ", ")),
                          indent = 2, exdent = 4))
    }
    if (other.pkgs && length(opkgver)) {
        if (is.null(pkgs))
            opkgver <- opkgver[sort(names(opkgver))]
        if (citations) {
            bibs <- write.bib(names(opkgver), file = file, append = TRUE,
                              verbose = FALSE)
            all_bibs <- c(all_bibs, bibs)
            key <- bibs$key
        }
        z <- c(z, "  \\item Used packages: ", "  \\begin{itemize}",
               formatPkgs(names(opkgver), opkgver, key), "  \\end{itemize}")
    }
    if (namespace.pkgs && length(nspkgver)) {
        nspkgver <- nspkgver[sort(names(nspkgver))]
        if (citations) {
            bibs <- write.bib(names(nspkgver), file = file, append = TRUE,
                              verbose = FALSE)
            all_bibs <- c(all_bibs, bibs)
            key <- bibs$key
        }
        z <- c(z, "  \\item Loaded via a namespace (and not attached): ",
               "  \\begin{itemize}",
               formatPkgs(names(nspkgver), nspkgver, key), "  \\end{itemize}")
    }
    z <- c(z, "\\end{itemize}")

    if (citations && !is.null(file)) {
        message("Written ", length(all_bibs), " BibTeX entries to file '", file,
                "' ...")
        message("Use \\bibliography{", file, "} to include citations.\n\n")
    }
    if (citations && is.null(file)) {
        attr(z, "BibTeX") <- all_bibs
        class(z) <- c("LatexBibtex", "Latex")
        return(z)
    } else {
        class(z) <- "Latex"
        return(z)
    }
}

print.LatexBibtex <- function(x, ...) {
    NextMethod("print", x)
    cat("\n\n")
    print(attr(x, "BibTeX"))
    invisible(x)
}

toLatex.LatexBibtex <- function(object, ...) {
    attributes(object) <- NULL
    class(object) <- "Latex"
    object
}

toBibtex.LatexBibtex <- function(object, ...) {
    object <- toBibtex(attr(object, "BibTeX"))
    object
}

formatPkgs <- function(name, vers, key, citecommand = "\\citep") {
    if (!is.null(key)) {
        key <- sapply(name, function(x)
            paste(key[grep(paste0("^pkg:", x, "[[:digit:]]*$"), key)],
                  collapse = ","))
        cites <- paste0(citecommand, "{", key, "}")
        cites[is.null(key)] <- ""
    } else {
        cites <- rep("", length(name))
    }
    paste0("\\item ", name, " (vers. ", vers, ") ", cites)
}

## modified version based on R package version 0.3-5.
##
## URL of original package: http://CRAN.R-project.org/package=bibtex
## Authors of R package bibtex (inlcuding write.bib):
##   Romain Francois, Kurt Hornik
## Licence of R package bibtex: GPL-2 | GPL-3
write.bib <- function(entry = "base", file = NULL,
                      append = FALSE, verbose = TRUE) {

    ## define bibs
    bibs <- if (inherits(entry, "bibentry")) {
        entry
    } else {
        if (length(entry) == 0) {
                if (verbose)
                    message("Empty package list: nothing to be done.")
                return(invisible(""))
        }
        if (is.character(entry)) {
            ## save names of packages
            pkgs <- entry
            bibs <- sapply(pkgs, function(x) citation(x), simplify = FALSE)
            n.installed <- length(bibs)
            ok <- sapply(bibs, inherits, "bibentry")
            pkgs <- pkgs[ok]
            bibs <- bibs[ok]
            n.converted <- sum(ok)
            ## generate unique keys
            pkgs <- lapply(seq_along(pkgs), function(i)
                           if (length(bibs[[i]]) > 1) {
                               paste0(pkgs[i], 1:length(bibs[[i]]))
                           } else {
                               pkgs[i]
                           })
            pkgs <- do.call("c", pkgs)
            bibs <- do.call("c", bibs)
            ## add keys to bibentries
            bibs <- mapply(function(b, k) {
                b$key <- paste0("pkg:", k)
                b
            }, bibs, pkgs, SIMPLIFY = FALSE)
            bibs <- do.call("c", bibs)
            if (verbose)
                message("Converted ", n.converted, " of ", n.installed,
                        " package citations to BibTeX")
            bibs
        } else {
            stop("Invalid argument 'entry': ",
                 "expected a bibentry object or a character vector ",
                 "of package names.")
        }
    }

    if (length(bibs) == 0) {
        if (verbose)
            message("Empty bibentry list: nothing to be done.")
        return(invisible())
    }
    if (!is.null(file)) {
        if (is.character(file)) {
            if (!grepl("\\.bib$", file))
                file <- paste(file, ".bib", sep = "")
        }
        fh <- file(file, open = ifelse(append, "a+", "w+"))
        on.exit(if (isOpen(fh)) close(fh))
        if (verbose)
            message("Writing ", length(bibs), " BibTeX entries ... ",
                    appendLF = FALSE)
        writeLines(toBibtex(bibs), fh)
        if (verbose)
            message("OK\nResults written to file '", file, "'")
        return(invisible(bibs))
    } else {
        return(bibs)
    }
}
