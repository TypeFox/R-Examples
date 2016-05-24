
labels2R <- function(lfile, rfile, dfname = "b", echo = FALSE)
{
    if (missing(lfile)) 
        stop(gettext("The name of file with labels is required.", 
            domain = "R-descr"))
    if (missing(rfile)) 
        stop(gettext("The name of file with R code to is required.", 
            domain = "R-descr"))
    if (!is.character(lfile)) 
        stop("lfile must be of class character.")
    if (!is.character(rfile)) 
        stop("rfile must be of class character.")
    if (!is.character(dfname)) 
        stop("dfname must be of class character.")
    infile <- path.expand(lfile[1])
    outfile <- path.expand(rfile[1])
    if (!file.exists(infile)) {
        msg <- paste(gettext("File not found:", domain = "R-descr"), 
            lfile)
        stop(msg)
    }
    if (file.exists(outfile)) {
        unlink(outfile)
    }
    input <- readLines(infile)
    nlines <- length(input)
    lnum <- 1
    while (lnum <= nlines) {
        cline <- input[lnum]
        if(echo)
            cat("[", lnum, "]", cline, "\n", sep = "")
        varname <- NULL
        varlab <- NULL
        lev <- NULL
        lab <- NULL
        exclud <- NULL
        if (cline != "" && grep("^[a-zA-Z]", cline) == 1) {
            varname <- sub("^([a-zA-Z0-9_\\.]*).*", "\\1", cline)
            if (grep(" ", cline) == 1) 
                varlab <- sub("(\\w|\\.|_)* (.*)", "\\2", cline)
            lnum <- lnum + 1
            cline <- input[lnum]
            if(echo)
                cat("[", lnum, "]", cline, "\n", sep = "")
            nlev = 0
            while (cline != "" && (grep("^[0-9]* ", cline) == 1 || grep("^-[0-9]* ", cline) == 1) && lnum <= nlines) {
                nlev <- nlev + 1
                lev[nlev] <- sub("^(-*[0-9]*) .*", "\\1", cline)
                lab[nlev] <- sub("^-*[0-9]* (.*)", "\\1", cline)
                if (lnum < nlines) {
                  lnum <- lnum + 1
                  cline <- input[lnum]
                  if(echo)
                      cat("[", lnum, "]", cline, "\n", sep = "")
                }
            }
            if (!is.null(lev)) 
                cat(dfname, "$", varname, " <- factor(", dfname, 
                  "$", varname, ",\n    levels = c(", paste(lev, 
                    collapse = ", "), "),\n    labels = c(\"", 
                  paste(lab, collapse = "\", \""), "\"))\n", 
                  sep = "", file = outfile, append = TRUE)
            if (!is.null(varlab)) 
                cat("attr(", dfname, "$", varname, ", \"label\") <- \"", 
                  varlab, "\"\n", sep = "", file = outfile, append = TRUE)
            cat("\n", file = outfile, append = TRUE)
        }
        if (lnum <= nlines)
            lnum <- lnum + 1
    }
    return(invisible(NULL))
}


