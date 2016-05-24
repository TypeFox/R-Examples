#' Function to generate a dose-response model from a data frame.
#'
#' @param x a data frame that has at least the columns listed in the \code{form}
#'   argument (e.g. "response" and "dose", see examples).
#' @param form a formula specifying the column names for the response and dose.
#'   Defaults to NULL.
#' @param model one of 4 options:
#' \itemize{
#'    \item LL.3 = Log Logistic 3 parameter with a lower limit of 0.
#'    \item LL.4 = Log Logistic 4 parameter with lower limit estimated.
#'    \item W1.4 = Weibul 4 parameter type 1.
#'    \item W2.4 = Weibul 4 parameter type 2.
#' }
#' @param idcol the name of the column that identifies the samples (case sensitive).
#'
#' @author Zhian N. Kamvar
#' @export
#' @keywords internal
#'
#' @details A wrapper function for \code{\link[drc]{drm}}, this will attempt to
#'   catch errors generated due to non-finite responses. When these are
#'   encounterd, a message that the model was not evaluated will be printed to
#'   the screen and NA will be returned.
get_drm <- function(x, form = NULL,
                    model = c("LL.3", "LL.4", "W1.4", "W2.4"), idcol = "ID"){
  ARGS <- c("LL.3", "LL.4", "W1.4", "W2.4")
  model <- match.arg(model, ARGS)
  if (model == "LL.3"){
    mod.names <- c("Slope", "Upper Limit", "ED50" )
  } else {
    mod.names <- c("Slope", "Lower Limit", "Upper Limit", "ED50" )
  }
  if (is.null(form)){
    the_call <- match.call()
    the_call[["form"]] <- response ~ dose
    the_call <- utils::capture.output(print(the_call))
    msg <- paste("please supply a formula.\n\nExample:\n\t", the_call)
    stop(msg)
  }
  MODEL <- match.fun(model)
  res <- tryCatch(drc::drm(form, data = x,
                           fct = MODEL(names = mod.names),
                           na.action = na.omit),
                  error = function(e) cat("Not evaluated:", x[[idcol]][1], "\n"),
                  finally = cat(""))
  return(res)
}
