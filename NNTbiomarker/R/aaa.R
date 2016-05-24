cat("======== aaa.R  ================\n")


#printFunctionBody = function(f) attributes(attributes(f)$srcref)$srcfile$lines

printFunctionBody = function(f)
  (deparse((body(f))))

### Convenience utilities borrowed from mvbutils

#' \%&\% string concatenation
#'
#' From mvbutils
#' @param a a string
#' @param b another string
#' @return paste0(a, b)
#'
`%&%` = function (a, b)
        paste(a, b, sep = "")

`%except%` = function (vector, condition)
	vector[match(vector, condition, 0) == 0]

"cq" = function (...)
{
    as.character(sapply(as.list(match.call(expand.dots = TRUE))[-1],
        as.character))
}

as.cat = function (x)
{
  stopifnot(is.character(x))
  oldClass(x) <- "cat"
  x
}

catn=function(...) cat(..., "\n")

withNames =
  function(x, n) {temp = data.frame(x=x,n=n);
                  x = temp$x;
                  n = temp$n;
                  names(x) <- n;
                  x}

# .onAttach = function(libname, pkgname){
#   cat("libname, pkgname:  ", libname, pkgname, "\n")
#   cat("search[1:2]: ", search()[1:2], "\n")
# }

clear = function(keep=c(".ctde", "startup", ".NNTbiomarker.verboseOptions")){
  answer <- repeat {
    cat("Delete ALL files in .GlobalEnv except ",
      paste(keep, collapse="&"), "?\n  (cannot be undone): ")
    answer <- readline()
    answer <- gsub("(\\w)", "\\U\\1", answer, perl=T)
    answer <- pmatch(answer, c("YES",  "NO", "N"))
      if (!is.na(answer)) {
        if(answer %in% 1)
        rm(list=setdiff(ls(all.names=T, pos=1), keep), pos=1)
      else
        cat("Aborted. No objects deleted.\n")
      return(invisible(NULL))
    }
  }
}


### other utilities
### inclusive , includes the endpoints
"%between%" = function(x, range) { (x<=range[2] & x>=range[1])}

