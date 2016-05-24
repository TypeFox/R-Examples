#' Get a list of implemented test functions with specific tags.
#'
#' @param tags [\code{character}]\cr
#'   Character vector of tags. All available tags can be determined with a call
#'   to \code{\link{getAvailableTags}}.
#' @param or [\code{logical(1)}]\cr
#'   Should all \code{tags} be assigned to the function or are single tags allowed
#'   as well?
#'   Default is \code{FALSE}.
#' @return [\code{character}]
#'   Named vector of function names with the given tags.
#' @examples
#' # list all functions which are unimodal
#' filterFunctionsByTags("unimodal")
#' # list all functions which are both unimodal and separable
#' filterFunctionsByTags(c("unimodal", "separable"))
#' # list all functions which are unimodal or separable
#' filterFunctionsByTags(c("multimodal", "separable"), or = TRUE)
#' @export
filterFunctionsByTags = function(tags, or = FALSE) {
  assertSubset(tags, choices = getAvailableTags(), empty.ok = FALSE)

  if (isSubset(c("single-objective", "multi-objective"), tags)) {
    stopf("Trying to search for both single- and multi-objective functions.")
  }

  fun.generators = getGeneratorFunctions()

  # helpers
  containsAllTags = function(tags, fun.tags) BBmisc::isSubset(tags, fun.tags)
  containsAtLeastOneTag = function(tags, fun.tags) any(tags %in% fun.tags)

  # what to check
  filter.fun = if (or) containsAtLeastOneTag else containsAllTags

  # filter by tags
  filtered.generators = Filter(function(fun) {
    fun.tags = attr(fun, "tags")
    return(filter.fun(tags, fun.tags))
  }, fun.generators)

  # some funs might occur multiple times in this case
  if (or) {
    filtered.generators = unique(filtered.generators)
  }

  # cleanup
  names = sapply(filtered.generators, function(fun) attr(fun, "name"))
  names(names) = NULL
  return(names)
}

# Get all generator objects.
#
# @return [function]
#   Vector of functions
getGeneratorFunctions = function() {
  # get all methods
  all.methods = ls("package:smoof")
  # get the function and not the names only
  all.methods = sapply(all.methods, get)
  # filter generators
  fun.generators = Filter(function(fun) inherits(fun, "smoof_generator"), all.methods)

  return(fun.generators)
}

getGeneratorByName = function(fun.name) {
  # get all methods of the package
  all.methods = sapply(ls("package:smoof"), get)
  fun.generator = Filter(function(fun) (fun.name %in% attr(fun, "name")), all.methods)
  if (length(fun.generator) == 0L) {
    return(NULL)
  }
  return(fun.generator[[1L]])
}
