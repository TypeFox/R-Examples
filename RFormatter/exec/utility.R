#!/usr/bin/Rscript

library(RFormatter)

main <- function() {
    # Format the R files given as command-line arguments.
    filenames <- commandArgs(TRUE)
    sapply(filenames, format_R_file)
}

format_R_file <- function(filename) {
    source <- paste(readLines(filename), collapse = "\n")
    formatR_arguments <- list(arrow = TRUE, width.cutoff = 80)
    formatted_source <- format_R_source_code(source, formatR_arguments)
    writeLines(formatted_source, filename)
}

if (!interactive()) {
    invisible(main())
}
