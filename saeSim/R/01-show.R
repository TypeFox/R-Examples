#' Show for \code{sim_setup}
#' 
#' This is the documentation for the show methods in the package \code{saeSim}. In case you don't know, \code{show } is for S4-classes like \code{print} for S3. If you don't know what that means, don't bother, there is no reason to call \code{show} directly, however there is the need to document it.
#' 
#' @inheritParams methods::show
#' 
#' @details Will print the head of a \code{sim_setup} to the console, after converting it to a \code{data.frame}.
#' @rdname showMethods
#' @export
setMethod("show", "sim_setup", function(object) {
  # this is essentially dplyr::print.tbl_df but there is no obvious 
  # constructor for tbl_df ...
  run <- sim_run_once(object)
  if (inherits(run, "data.frame")) {
    dat <- as.data.frame(run)
    cat("data.frame ", dim_desc(dat), "\n", sep = "")
    cat("\n")
    print(trunc_mat(dat, n = 6, width = NULL))
    invisible(dat)
  } else if (isS4(run)) {
    show(run)
    invisible(run)
  } else {
    print(run)
  }
})

