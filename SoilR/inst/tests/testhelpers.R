#!/usr/bin/Rscript
# vim:set ff=unix expandtab ts=2 sw=2:
checkWarning <- function(expr,silent=TRUE) {
    checkTrue(
    inherits(
       tryCatch(
         eval(expr, envir = parent.frame()), 
         silent = silent,
         warning=function(w){w}
       ),
       "warning"
    ),
        "Warning not generated as expected\n")
}             
