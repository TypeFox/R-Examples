#' Short \code{selection} summary
#'
#' \code{\link{selection}} summary
#' @param x \code{selection} object.
#' @param \ldots Other options.
#' @return The function returns the best subset of size q and its information
#'   criterion value. In the case of \code{seconds=TRUE} this information is
#'   returned for each alternative model.
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' @seealso \code{\link{selection}}.
#' @examples
#' library(FWDselect)
#' data(diabetes)
#' x = diabetes[ ,2:11]
#' y = diabetes[ ,1]
#' obj1 = selection(x, y, q = 1, method = "lm", criterion = "variance", cluster = FALSE)
#' obj1
#' @export



print.selection <- function(x = model, ...) {

  if (inherits(x, "selection")) {

    model <- x

    cat("\n")
    cat("****************************************************")

    cat("\nBest subset of size q =", length(model$Variable_numbers),
        ": ")
    cat(format(model$Variable_names))
    cat("\n")

    cat("\nInformation Criterion Value -", model$ic,
        ": ")
    cat(format(model$Information_Criterion))
    cat("\n")

    cat("****************************************************")
    cat("\n")
    if (model$seconds == T) {

      cont = 0
      auxic <- c()
      auxmodel <- c()
      for (i in 1:model$nmodels) {

        if (i == 1) {
          cont = cont + 10
        } else {
          cont = cont + 5
        }

        # change order (due to the cv)
        auxic[i] <- model[cont + 2]
        auxmodel[i] <- model[cont]
      }

      ii <- order(unlist(auxic))

      for (i in 1:model$nmodels) {

        cat("\nAternative (", i, ") subset of size q =",
            length(model$Variable_numbers),
            ": ")
        cat(format(auxmodel[ii[i]]))

        cat("\n")
        cat("\nInformation Criterion Value -",
            model$ic, ": ")
        cat(format(auxic[ii[i]]))
        cat("\n")
        cat("\n")
      }
    }
  }else{
    stop("Argument x must be either selection object.")
  }
}







