
## Method to overwrite parameters of an existing MaxControl object
addControlDddot <- function(x, ...) {
   ## add ... to the control
   dddot <- list(...)
   dddot <- dddot[names(dddot) %in% openParam(x)]
   addControlList(x, dddot)
   ##
}

setMethod("maxControl", "MaxControl", addControlDddot)
