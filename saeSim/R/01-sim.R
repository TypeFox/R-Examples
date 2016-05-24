#' Start simulation
#' 
#' This function will start the simulation. Use the printing method as long as
#' you are testing the scenario.
#' 
#' @param x a \code{sim_setup}
#' @param R number of repetitions.
#' @param path optional path in which the simulation results can be saved. They
#'   will we coerced to a \code{data.frame} and then saved as 'csv'.
#' @param overwrite \code{TRUE}/\code{FALSE}. If \code{TRUE} files in
#'   \code{path} are replaced. If \code{FALSE} files in \code{path} are not
#'   replaced and simulation will not be recomputed.
#' @param ... arguments passed to \code{\link{parallelStart}}.
#' @param libs arguments passed to \code{\link{parallelLibrary}}. Will be used
#'   in a call to \code{\link{do.call}} after coersion with
#'   \code{\link{as.list}}.
#' @param exports arguments passed to \code{\link{parallelExport}}. Will be used
#'   in a call to \code{\link{do.call}} after coersion with
#'   \code{\link{as.list}}.
#' @param suffix an optional suffix of file names.
#' @param fileExt the file extension. Default is ".csv" - alternative it can be
#'   ".RData".
#' 
#' @details The package parallelMap is utilized as back-end for parallel computations.
#' 
#' Use the argument \code{path} to store the simulation results in a directory.
#' This may be a good idea for long running simulations and for those using
#' large \code{data.frame}s. You can use \code{\link{sim_read_data}} to read
#' them in. The return value will change to NULL in each run.
#'
#' @return The return value is a list. The elements are the results of each
#'   simulation run, typically of class \code{data.frame}. In case you specified
#'   \code{path}, each element is \code{NULL}.
#' 
#' @rdname sim
#' @export
#' @examples
#' setup <- sim_base_lm()
#' resultList <- sim(setup, R = 1)
#' 
#' # For parallel computations you may need to export objects
#' localFun <- function() cat("Hello World!")
#' comp_fun <- function(dat) {
#'   localFun()
#'   dat
#' }
#' 
#' res <- sim_base_lm() %>% 
#'   sim_comp_pop(comp_fun) %>% 
#'   sim(R = 2, 
#'       mode = "socket", cpus = 2,
#'       exports = "localFun")
#' 
#' str(res)
sim <- function(x, R = 1, path = NULL, overwrite = TRUE, ..., suffix = NULL, fileExt = ".csv", libs = NULL, exports = NULL) {
  
  parallelStart(...)
  do.call(parallelLibrary, as.list(libs))
  do.call(parallelExport, as.list(exports))
  res <- parallelLapply(
    1:R, map_fun, object = x, 
    path = path, overwrite = overwrite, suffix = suffix, fileExt = fileExt
  )
  parallelStop()
  res
  
}

map_fun <- function(i, object, path, overwrite, suffix, fileExt) {
  filename <- make_sim_filename(i, object, path, suffix, fileExt)
  res <- if (needs_recompute(filename, overwrite)) {
    sim_run_it(object, i)
  } else {
    NULL
  }
  sim_write_results(res, path, filename, fileExt)
}

make_sim_filename <- function(i, object, path, suffix, fileExt) {
  suffix <- if (is.null(suffix)) "" else paste0("-", suffix)
  if (!is.null(path)) paste0(path, "/", object@simName, i, suffix, fileExt) else NULL
}

needs_recompute <- function(filename, overwrite) {
  if (is.null(filename)) return(TRUE) # path = NULL
  if (!file.exists(filename)) return(file.create(filename)) # always compute if it doesn't exist 
  else return(overwrite)
}

sim_run_it <- function(object, i) {
  res <- sim_run_once(object)
  res$idR <- i
  res$simName <- object@simName
  res
}

sim_run_once <- function(x) {
  Reduce(function(x, f) f(x), x, x@base)
}

sim_write_results <- function(res, path, filename, fileExt) {
  if (!is.null(path) && !is.null(res)) {
    if (identical(fileExt, ".csv")) {
      res <- as.data.frame(res)
      write.csv(res, file = filename, row.names = FALSE)
    } else if (identical(fileExt, ".RData")) {
      save(list = "res", file = filename)
    } else {
      stop("Unknown file extension for this result.")
    }
    NULL
  } else {
    res
  }
}

