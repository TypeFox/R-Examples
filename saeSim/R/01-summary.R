#' Summary for a sim_setup
#' 
#' Reports a summary of the simulation setup.
#' 
#' @param object a \code{sim_setup}.
#' @param ... has no effect.
#' 
#' @export
#' 
#' @examples
#' summary(sim_base_lm())
setMethod("summary", c(object = "sim_setup"), function(object, ...) {
  callList <- lapply(S3Part(object, strictS3=TRUE), getCalls)
  expr <- parse(text = do.call(paste, c(callList, sep = " %>%\n\t\t")))
  time <- system.time(dat <- as.data.frame(object))
  dimension <- dim(dat)
  
  new("summary.sim_setup", 
      sim_setup = object,
      dim = dimension, 
      duration = as.table(time), 
      expression = expr)
})

getCalls <- function(simFun) {
  cl <- slot(simFun, "call")
  cl$simSetup <- NULL
  paste(deparse(cl), collapse = "\n")
} 

#' @rdname showMethods
#' @inheritParams methods::show
#' @export
setMethod("show", "summary.sim_setup", function(object) {
  cat("General Information about", object@sim_setup@simName, "simulation set-up:\n")
  print(object@expression)
  
  tmp <- object@duration
  cat("\nApproximating the expected duration:\n")
  cat("A single run takes ")
  cat(tmp["elapsed"], "seconds. 500 * ", tmp["elapsed"], "=", 
      500 * tmp["elapsed"], "seconds.\n")
})



