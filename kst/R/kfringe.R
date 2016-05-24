### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kfringe.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets), library(relations)
###
### 2009-08-06: created
###

kfringe <- function(x, state=NULL, operation=c("inner", "outer")) {

   ### check type
   operation <- match.arg(operation)

   ### check x
   if (!inherits(x, "kstructure"))
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))


   ### check state
   if (!is.null(state) && !set_contains_element(x, state))
      stop(sprintf("Specified state is no element of %s", dQuote("x")))

   ### return results
   if (operation=="inner")
      kfringe_inner(x, state=state)
   else
      kfringe_outer(x, state=state)


}
