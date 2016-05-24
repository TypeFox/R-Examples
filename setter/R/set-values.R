#' Set attributes of variables
#'
#' Set attributes of a variable, then return that variable.
#' @param x An R variable.
#' @param value A value to set the attribute to.
#' @param from A variable to take the attribute value from.
#' @param ... Name-value pairs of attributes to set.
#' @param .dots A named list of attributes to set.
#' @param attribs A character vector of attributes to copy.
#' @return \code{x} is returned, with updated attributes.
#' @seealso \code{\link[base]{attributes}}, \code{\link[base]{attr}}
#' @name setter
#' @aliases setters mutator
#' @examples
#' if(requireNamespace("magrittr"))
#' {
#'   `%>%` <- magrittr::`%>%`
#'   # Convert a vector to a matrix by setting the dimensions and their names.
#'   m <- 1:12 %>%
#'     set_dim(3:4) %>%
#'     set_dimnames(list(letters[1:3], LETTERS[1:4])) %>%
#'     print
#'
#'   # Copy attributes from one variable to another using copy_* fns.
#'   month.abb %>%
#'     copy_dim(m) %>%
#'     copy_dimnames(m) %>%
#'     print
#'
#'   # Same again, using copy_attributes
#'   month.abb %>%
#'     copy_attributes(m, c("dim", "dimnames"))
#'
#'   # Same again, in this case you can copy most/all the attributes from m.
#'   month.abb %>%
#'     copy_most_attributes(m)
#'   month.abb %>%
#'     copy_all_attributes(m)
#'
#'   # To quickly convert a list into a data.frame, set the class and row names.
#'   list(a = (1:5) ^ 2, b = pi ^ (1:5)) %>%
#'     set_class("data.frame") %>%
#'     set_rownames() %>%   # data.frames have a default
#'     print
#'
#'   # Or equivalently, using attributes
#'   list(a = (1:5) ^ 2, b = pi ^ (1:5)) %>%
#'     set_attributes(class = "data.frame", row.names = .set_row_names(5)) %>%
#'     print
#' } else
#' {
#'   message('This example requires the magrittr package.  Please run install.packages("magrittr").')
#' }
#' @importFrom assertive.base coerce_to
NULL

# Class & mode ------------------------------------------------------------

#' @rdname setter
#' @export
set_class <- function(x, value)
{
  value <- coerce_to(value, "character")
  class(x) <- value
  x
}

#' @rdname setter
#' @export
copy_class <- function(x, from)
{
  class(x) <- class(from)
  x
}

#' @rdname setter
#' @export
set_mode <- function(x, value)
{
  value <- coerce_to(value, "character")
  mode(x) <- value
  x
}

#' @rdname setter
#' @export
copy_mode <- function(x, from)
{
  mode(x) <- mode(from)
  x
}

#' @rdname setter
#' @export
set_storage_mode <- function(x, value)
{
  value <- coerce_to(value, "character")
  storage.mode(x) <- value
  x
}

#' @rdname setter
#' @export
copy_storage_mode <- function(x, from)
{
  storage.mode(x) <- storage.mode(from)
  x
}

# Dimensions & length -----------------------------------------------------

#' @rdname setter
#' @export
set_dim <- function(x, value)
{
  value <- coerce_to(value, "integer")
  dim(x) <- value
  x
}

#' @rdname setter
#' @export
copy_dim <- function(x, from)
{
  dim(x) <- dim(from)
  x
}

#' @rdname setter
#' @export
set_length <- function(x, value)
{
  value <- coerce_to(value, "integer")
  length(x) <- value
  x
}

#' @rdname setter
#' @export
copy_length <- function(x, from)
{
  length(x) <- length(from)
  x
}

# Names -------------------------------------------------------------------

#' @rdname setter
#' @export
set_names <- function(x, value)
{
  value <- coerce_to(value, "character")
  names(x) <- value
  x
}

#' @rdname setter
#' @export
copy_names <- function(x, from)
{
  names(x) <- names(from)
  x
}

#' @rdname setter
#' @export
set_colnames <- function(x, value)
{
  value <- coerce_to(value, "character")
  colnames(x) <- value
  x
}

#' @rdname setter
#' @export
copy_colnames <- function(x, from)
{
  colnames(x) <- colnames(from)
  x
}

#' @rdname setter
#' @export
set_rownames <- function(x, value)
{
  # data.frames can use magic from .set_row_names to provide a default
  # data.tables and tibble::data_frames will mostly ignore row names.
  x_is_df <- is.data.frame(x)
  if(x_is_df && missing(value) || is.null(value))
  {
    # dim and nrow may lie if you've just list() %>% set_class("data.frame")
    attr(x, "row.names") <- .set_row_names(length(x[[1]]))
  } else if(x_is_df && is.numeric(value) && length(value) == 1)
  {
    attr(x, "row.names") <- .set_row_names(value)
  } else
  {
    value <- coerce_to(value, "character")
    rownames(x) <- value
  }
  x
}

#' @rdname setter
#' @export
copy_rownames <- function(x, from)
{
  rownames(x) <- rownames(from)
  x
}

#' @rdname setter
#' @export
set_dimnames <- function(x, value)
{
  # dimnames<- should usually take a list of character vectors, but it
  # has its own complicated set of rules for dealing with weird input
  # so don't check here; just let it do its own thing.
  dimnames(x) <- value
  x
}

#' @rdname setter
#' @export
copy_dimnames <- function(x, from)
{
  dimnames(x) <- dimnames(from)
  x
}

# Levels ------------------------------------------------------------------

#' @rdname setter
#' @export
set_levels <- function(x, value)
{
  value <- coerce_to(value, "character")
  levels(x) <- value
  x
}

#' @rdname setter
#' @export
copy_levels <- function(x, from)
{
  levels(x) <- levels(from)
  x
}

# Comments ----------------------------------------------------------------

#' @rdname setter
#' @export
set_comment <- function(x, value)
{
  value <- coerce_to(value, "character")
  comment(x) <- value
  x
}

#' @rdname setter
#' @export
copy_comment <- function(x, from)
{
  comment(x) <- comment(from)
  x
}

# Parts of functions ------------------------------------------------------
# Not working with pipes; exclude for now.

#  @param envir The environment in which the function should be defined.  See
#  \code{\link[base]{formals}} or \code{\link[base]{body}}.

# #' @rdname setter
# #' @export
# set_formals <- function(x, value, envir = environment(x))
# {
#   # formals<- has its own complicated handling of inputs; don't try to replicate
#   force(envir)
#   formals(x, envir) <- value
#   x
# }
#
# #' @rdname setter
# #' @export
# copy_formals <- function(x, from, envir = environment(x))
# {
#   force(envir)
#   formals(x, envir) <- formals(from)
#   x
# }
#
# #' @rdname setter
# #' @export
# set_body <- function(x, value, envir = environment(x))
# {
#   force(envir)
#   body(x, envir) <- value
#   x
# }

# #' @rdname setter
# #' @export
# copy_body <- function(x, from, envir = environment(x))
# {
#   force(envir)
#   body(x, envir) <- body(from)
#   x
# }



# General attributes ------------------------------------------------------

#' @importFrom assertive.base merge_dots_with_list
#' @rdname setter
#' @export
set_attributes <- function(x, ..., .dots = list())
{
  values <- merge_dots_with_list(..., l = .dots)
  attributes <- names(values)
  for(i in seq_along(attributes))
  {
    attr(x, attributes[i]) <- values[[i]]
  }
  x
}

#' @rdname setter
#' @importFrom assertive.base coerce_to
#' @export
copy_attributes <- function(x, from, attribs)
{
  attribs <- coerce_to(attribs, "character")
  for(attribute in attribs)
  {
    attr(x, attribute) <- attr(from, attribute)
  }
  x
}

#' @rdname setter
#' @export
copy_all_attributes <- function(x, from)
{
  attributes(x) <- attributes(from)
  x
}

#' @rdname setter
#' @export
copy_most_attributes <- function(x, from)
{
  mostattributes(x) <- attributes(from)
  x
}
