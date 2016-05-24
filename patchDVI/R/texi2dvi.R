# This is extracted from the file in the original copyright notice below,
# and modified for patchDVI.

#  File src/library/tools/R/utils.R
#  Part of the R package, http://www.R-project.org
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


## x without y, as in the examples of ?match.
`%w/o%` <-
function(x, y)
    x[!x %in% y]

Rtexinputs <- function() {
    if (getRversion() < "2.12.0")
    	file.path(R.home("share"), "texmf")
    else
	file.path(R.home("share"), "texmf", "tex", "latex")
}

texi2dvi <-
function(file, pdf = FALSE, clean = FALSE, quiet = FALSE,
         texi2dvi = getOption("texi2dvi"),
         texinputs = NULL, index = TRUE,
         links = NULL )
{
    ## Run texi2dvi on a latex file, or emulate it.

    if(is.null(texi2dvi) || !nzchar(texi2dvi))
        texi2dvi <- Sys.which("texi2dvi")

    envSep <- .Platform$path.sep
    texinputs0 <- texinputs
    
    Rtexmf <- file.path(R.home("share"), "texmf")
    
    ## "" forces use of default paths.
    texinputs <- paste(c(texinputs0, Rtexinputs(), ""), 
                      collapse = envSep)
    ## not clear if this is needed, but works
    if(.Platform$OS.type == "windows")
        texinputs <- gsub("\\", "/", texinputs, fixed = TRUE)
    if (getRversion() < "2.12.0")
        Rbstinputs <- Rtexmf
    else 
    	Rbstinputs <- file.path(Rtexmf, "bibtex", "bst")
    bstinputs <- paste(c(texinputs0, Rbstinputs, ""), 
                      collapse = envSep)

    otexinputs <- Sys.getenv("TEXINPUTS", unset = NA)
    if(is.na(otexinputs)) {
        on.exit(Sys.unsetenv("TEXINPUTS"))
        otexinputs <- "."
    } else on.exit(Sys.setenv(TEXINPUTS = otexinputs))
    Sys.setenv(TEXINPUTS = paste(otexinputs, texinputs, "", sep = envSep))
    bibinputs <- Sys.getenv("BIBINPUTS", unset = NA)
    if(is.na(bibinputs)) {
        on.exit(Sys.unsetenv("BIBINPUTS"), add = TRUE)
        bibinputs <- "."
    } else on.exit(Sys.setenv(BIBINPUTS = bibinputs, add = TRUE))
    Sys.setenv(BIBINPUTS = paste(bibinputs, texinputs, sep = envSep))
    obstinputs <- Sys.getenv("BSTINPUTS", unset = NA)
    if(is.na(obstinputs)) {
        on.exit(Sys.unsetenv("BSTINPUTS"), add = TRUE)
        obstinputs <- "."
    } else on.exit(Sys.setenv(BSTINPUTS = obstinputs), add = TRUE)
    Sys.setenv(BSTINPUTS = paste(obstinputs, bstinputs, sep = envSep))

    if(index && nzchar(texi2dvi) && .Platform$OS.type == "windows") { # MiKTeX on Windows
        extra <- ""
        if (is.null(links))
            opt_links <- if(pdf) "--tex-option=-synctex=-1" else "--tex-option=--src-specials"
        else
            opt_links <- links
        ext <- if(pdf) "pdf" else "dvi"
        opt_pdf <- if(pdf) "--pdf" else ""
        file.create(".timestamp")
        opt_quiet <- if(quiet) "--quiet" else ""

        ## look for MiKTeX (which this almost certainly is)
        ## and set the path to R's style files.
        ## -I works in MiKTeX >= 2.4, at least
        ## http://docs.miktex.org/manual/texify.html
        cmd <- paste(shQuote(texi2dvi), "--version")
        if (!quiet) message(cmd, "\n")
        ver <- system(cmd, intern = TRUE)
        if(length(grep("MiKTeX", ver[1L]))) {
            ## AFAICS need separate -I for each element of texinputs.
            texinputs <- c(texinputs0, 
                           Rtexinputs(),
                           if (getRversion() < "2.12.0") NULL else Rbstinputs)
            paths <- paste ("-I", shQuote(texinputs))
            extra <- paste(extra, paste(paths, collapse = " "))
        }
        ## this only gives a failure in some cases, e.g. not for bibtex errors.
        cmd <- paste(shQuote(texi2dvi), opt_quiet, opt_pdf, opt_links,
                     shQuote(file), extra)
        if (!quiet) message(cmd, "\n")
        consoleLog <- system(cmd, intern=TRUE, ignore.stderr=TRUE)
        if(clean) {
            out_file <- paste(tools::file_path_sans_ext(file), ext, sep = ".")
            files <- list.files(all.files = TRUE) %w/o% c(".", "..",
                                                          out_file)
            file.remove(files[file_test("-nt", files, ".timestamp")])
        }
        file.remove(".timestamp")
    } else {
        ## Do not have Synctex-compatible texi2dvi or don't want to index
        ## Needed everywhere except for MiKTeX

        ## If it is called with MiKTeX then TEXINPUTS etc will be ignored.

        if (is.null(links))
            opt_links <- if(pdf) "--synctex=-1" else "--src-specials"
        else
            opt_links <- links
        texfile <- shQuote(file)
        base <- tools::file_path_sans_ext(file)
        idxfile <- paste(base, ".idx", sep="")
        latex <- if(pdf) Sys.getenv("PDFLATEX", "pdflatex")
        else  Sys.getenv("LATEX", "latex")
        bibtex <- Sys.getenv("BIBTEX", "bibtex")
        makeindex <- Sys.getenv("MAKEINDEX", "makeindex")
        cmd <- paste(shQuote(latex), "-interaction=nonstopmode", opt_links, texfile)
        if (!quiet) message(cmd, "\n")
        consoleLog <- system(cmd, intern = TRUE)
        status <- attr(consoleLog, "status")
        if(!is.null(status) && status)
            warning(gettextf("unable to run '%s' on '%s'", latex, file),
                 domain = NA)
        nmiss <- length(grep("^LaTeX Warning:.*Citation.*undefined",
                             readLines(paste(base, ".log", sep = ""))))
        for(iter in 1L:10L) { ## safety check
            ## This might fail as the citations have been included in the Rnw
            if(nmiss) system(paste(shQuote(bibtex), shQuote(base)))
            nmiss_prev <- nmiss
            if(index && file.exists(idxfile)) {
                cmd <- paste(shQuote(makeindex), shQuote(idxfile))
                if (!quiet) message(cmd, "\n")
                if(system(cmd))
                    stop(gettextf("unable to run '%s' on '%s'",
                                  makeindex, idxfile),
                         domain = NA)
            }
            cmd <- paste(shQuote(latex), "-interaction=nonstopmode", opt_links, texfile)
            if (!quiet) message(cmd, "\n")
	    consoleLog <- system(cmd, intern = TRUE)
	    status <- attr(consoleLog, "status")
            if(!is.null(status) && status)
                warning(gettextf("unable to run %s on '%s'", latex, file), domain = NA)
            Log <- readLines(paste(base, ".log", sep = ""))
            nmiss <- length(grep("^LaTeX Warning:.*Citation.*undefined", Log))
            if(nmiss == nmiss_prev &&
               !length(grep("Rerun to get", Log)) ) break
        }
    }
    invisible(consoleLog)
}

