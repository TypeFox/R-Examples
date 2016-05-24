#' Creates a Makefile rule
#'
#' A rule in a \code{Makefile} consists of a (list of) targets which may
#' depend on one or more dependencies each.  Optionally, a script is executed to
#' create the target.  Generally, multiple targets mean that the rule is
#' identical for each of the individual targets, and multiple dependencies mean
#' that \emph{all} of them are required to build \emph{each} of the targets.
#' In the script, the target can be referenced by \code{$@@}, and the first
#' dependency can be referenced by \code{$<}.  Note that the dollar sign has a
#' special meaning in a \code{Makefile}, use \code{$$} in scripts that need
#' to use the dollar sign themselves.
#'
#' Use the
#' \code{\link[base]{c}} function or the \code{\link[base]{+}} operator
#' to append rules to groups and Makefiles.
#'
#' @param targets \code{[character]}\cr Target names
#' @param deps \code{[character]}\cr Dependency names
#' @param script \code{[character]}\cr A script to execute to build the targets.
#' @return An object of class \code{MakefileR_rule}
#' @seealso \code{\link{makefile}}
#' @family items
#'
#' @examples
#' make_rule("all", c("first_target", "second_target"))
#' make_rule(".FORCE")
#' make_rule("first_target", ".FORCE", "echo 'Building first target'")
#' make_rule("second_target", "first_target",
#'  c("echo 'Building second target'", "echo 'Done'"))
#'
#' makefile() +
#'   make_rule("all", c("first_target", "second_target")) +
#'   make_rule(".FORCE") +
#'   make_rule("first_target", ".FORCE", "echo 'Building first target'") +
#'   make_rule("second_target", "first_target",
#'     c("echo 'Building second target'", "echo 'Done'"))
#'
#' @references \url{https://www.gnu.org/software/make/manual/}
#'
#' @export
make_rule <- function(targets, deps = NULL, script = NULL) {
  if (length(targets) == 0L) {
    stop("At least one target is required.")
  }
  structure(
    list(
      targets = targets,
      deps = deps,
      script = script
    ),
    class = c("MakefileR_rule", "MakefileR"))
}

#' @export
format.MakefileR_rule <- function(x, ...) {
  c(
    sprintf("%s:%s",
            combine_targets(x$targets),
            if (!is.null(x$deps))
              combine_targets(c("", x$deps))
            else
              ""),
    if (!is.null(x$script)) paste0("\t", x$script)
  )
}

combine_targets <- function(targets) {
  paste(targets, collapse = " ")
}
