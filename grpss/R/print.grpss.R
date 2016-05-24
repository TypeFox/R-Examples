#' Print the \code{grpss} object
#' @description Prints the "\code{grpss}" object obtained from \code{\link{grpss}}.
#' @param x A fitted "\code{grpss}" object.
#' @param ... Not used.
#' @details This function prints out the group screening results obtained from
#' \code{\link{grpss}} when "\code{selection = FALSE}".
#' @return The \code{call}, \code{criterion}, \code{threshold} and \code{screened groups} will
#' be printed out.
#' @author Debin Qiu, Jeongyoun Ahn
#' @seealso \code{\link{grpss}}
#' @export
print.grpss <- function(x,...) {
  cat("Call: \n")
  print(x$call)
  criteria <- switch(x$criterion, gSIS = "SIS", gHOLP = "HOLP",
                     gAR2 = "AR2", gDC = "DC")
  cat("\nCriterion: group", criteria)
  cat("\nThreshold (ncut): ", x$ncut)
  cat("\nScreened groups: ", x$group.screen)
}
