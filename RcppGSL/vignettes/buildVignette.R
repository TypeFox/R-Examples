#!/usr/bin/Rscript

## use given argument(s) as target files, or else default to .Rnw files in directory
files <- if (length(commandArgs(TRUE)) == 0) dir(pattern="*.Rnw") else commandArgs(TRUE)

## convert all files from Rnw to pdf using the highlight driver
invisible(sapply(files, function(srcfile) {
    Sweave(srcfile, driver=highlight::HighlightWeaveLatex(boxes=TRUE))
    tools::texi2pdf(gsub(".Rnw", ".tex", srcfile))
}))
