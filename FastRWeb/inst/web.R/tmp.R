## this allows generated scripts to pass references to
## temporary files created in the scripts as URLs
run <- function(file, mime="text/html",...)
    WebResult("tmpfile", gsub("/", ".", file, fixed=TRUE), mime)
