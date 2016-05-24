#' Convert network to data frame.
#'
#' @note If the instance contains of \eqn{n} depots, the depot coordinates fill the
#'   first \eqn{n} rows of the data frame.
#' @template arg_network
#' @param row.names [\code{character}]\cr
#'   Row names for the result. Default is \code{NULL}.
#' @param optional [any]\cr
#'   Currently not used.
#' @param include.extras [\code{logical(1)}]\cr
#'   Include additional information like cluster membership and node type as specific columns?
#'   Default is \code{TRUE}.
#' @param ... [any]\cr
#'   Currently not used.
#' @return [\code{data.frame}]
#' @export
as.data.frame.Network = function(x,
  row.names = NULL,
  optional = FALSE,
  include.extras = TRUE,
  ...) {
  n = getNumberOfNodes(x)
  res = as.data.frame(x$coordinates)

  # depot coordinates always the first ones
  if (hasDepots(x)) {
    res = rbind(as.data.frame(x$depot.coordinates), res)
  }
  colnames(res) = paste0("x", seq(ncol(res)))

  assertFlag(include.extras)
  if (!is.null(row.names)) {
    assertCharacter(row.names, len = n, any.missing = FALSE)
  }

  if (include.extras) {
    if (!hasDepots(x)) {
      res$types = "node"
    } else {
      res$types = c(rep("depot", getNumberOfDepots(x)), rep("node", n))
    }
    if (!is.null(x$membership)) {
      if (!hasDepots(x)) {
        res$membership = x$membership
      } else {
        res$membership = c(rep(0, getNumberOfDepots(x)), x$membership)
      }
    }
  }
  as.data.frame(res, row.names = row.names, optional = optional, ...)
}
