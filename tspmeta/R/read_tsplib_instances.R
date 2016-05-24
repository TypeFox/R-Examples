#' Read in multiple TSPLIB style Traveling Salesman Problems
#' from a directory.
#'
#' @param path [\code{character(1)}]\cr
#'   Character string containing path to file in TSPLIB format.
#' @param pattern [\code{character(1)}]\cr
#'   Pattern of files under \code{path} that are considered as instances.
#' @param max_size [\code{numeric(1)}]\cr
#'   Upper bound for instance size (i.e. number of cities). Only applicable,
#'   if instance size is contained in file name. Default value ist 1000.
#' @param use_names [\code{logical(1)}]\cr
#'   Use base names of files as names of instances in returned list.
#' @param on_no_coords [\code{character(1)}]\cr
#'   How to handle instances which do not have any
#'   coordinates. Possible values are, \dQuote{stop} and \dQuote{warn}
#'   which either stop or raise a warning respectivly.
#'
#' @return A \code{list} List of \code{tsp_instance} objects.
#' @export
read_tsplib_instances = function(path, pattern = "*.tsp",
  max_size = 1000, use_names = TRUE, on_no_coords = "stop") {

  assertCharacter(path, len = 1L, any.missing = FALSE)
  assertCharacter(pattern, len = 1L, any.missing = FALSE)
  max_size = convertInteger(max_size)
  assertInteger(max_size, len = 1L, any.missing = FALSE)
  assertFlag(use_names, na.ok = FALSE)
  assertChoice(on_no_coords, choices = c("stop", "warn"))

  file_names = list.files(path, pattern = pattern, full.names = TRUE)
  instances = list()
  errs = data.frame(character(0), character(0), character(0))
  for (file_name in file_names) {
    bn = basename(file_name)
    #FIXME: errorhandling if this goes wrong!
    k = as.integer(str_extract(bn, "\\d+"))
    if (k <= max_size) {
      x = try(read_tsplib_instance(file_name))
      if (is.error(x)) {
        errs = rbind(errs, data.frame(bn, file_name, as.character(x)))
        warningf("Instance could not be parsed correctly: %s", bn)
      }	else {
        instances[[length(instances) + 1]] = x
        #FIXME: should be indexed with length(instances) + 1?
        if (use_names)
          names(instances)[length(instances)] = bn
      }
    }
  }

  list(
    instances = instances,
    errors = setColNames(errs, c("filename", "path", "msg"))
  )
}
