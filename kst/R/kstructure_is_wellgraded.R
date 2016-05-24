### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kstructure_is_wellgraded.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2009-05-06: created
###

kstructure_is_wellgraded <- function(x) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### check for well-gradedness
   siw <- TRUE
   lp <- lpath_is_gradation(lpath(x))
   if (!all(unlist(lp))) {
      siw <- FALSE
   } else {
      ifringe <- kfringe_inner(x, state=NULL)
      ofringe <- kfringe_outer(x, state=NULL)
      fringe <- list()
      for (i in seq_along(ifringe)) {
         fringe[[i]] <- tuple(ifringe[[i]], ofringe[[i]])
      }
      for (j in seq_len(length(fringe)-1)) {
         fri <- fringe[(j+1):length(fringe)]
         if (any(sapply(fri,function(z)all(z==fringe[[j]]), simplify=TRUE))) {
            siw <- FALSE
            break
         }
      }
   }

   ### return result
   siw

}
