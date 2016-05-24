#' Retrieve a data.frame result of univariate regression
#'
#' Return a concentrated result of univariate regression models.

#' @param x A reg object
#' @param save A logical, whether to save the concentrated result
#' @param file A character string naming a file, see \code{\link{write.table}} for more details, the default filepath is current working directory and the filename is the current time.
#' @param sep,row.names,\dots See \code{\link{write.table}} for more details
#' @importFrom utils write.table
#' @export
#' @seealso \code{\link{reg}}
#' @examples
#' reg_glm<-reg(data = diabetes, x=c(1:4),y = 5,
#' factor = c(1, 3, 4), model = 'glm')
#'
#' dataframe(reg_glm)
#' # dataframe(reg_glm, save = TRUE)
#' # dataframe(reg_glm, file = "C:/reg_glm_out.csv")

dataframe.reg <- function(x, save = FALSE, file = NULL, sep = ",", row.names = FALSE,
    ...) {
    if (class(x) != "reg") {
        stop("x should be a `reg` object.", call. = FALSE)
    }
    if (is.null(file)) {
        time <- format(Sys.time(), "%Y-%m-%d %H-%M")
        file_new <- paste0("reg_", time, ".csv")
    }
    result <- x$dataframe

    if (save) {
      if (!is.null(file))  {
        write.table(result, file = file, sep = sep, row.names = row.names,
                    ...)
        } else write.table(result, file = file_new, sep = sep, row.names = row.names,
                           ...)

    }
    return(result)
}

