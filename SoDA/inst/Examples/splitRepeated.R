splitRepeated <- function(lines,
           separator = "[[:space:]]+",
           numeric = NA) {
    value <- strsplit(lines, separator)
    if(!identical(numeric, FALSE)) {
       warned <- FALSE
       opt <- options(warn = -1); on.exit(options(opt))
       nValue <-  withCallingHandlers(lapply(value, as.numeric),
                warning = function(cond) warned <<- TRUE)
       if(!warned || identical(numeric, TRUE))
           value <- nValue
       options(opt)
       if(warned && identical(numeric, TRUE))
         warning("NAs introduced in coercing to numeric")
   }
    value
}
