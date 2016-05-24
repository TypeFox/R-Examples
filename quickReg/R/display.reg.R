#' @title summary the result of univariate regression
#'
#' @description return summary information of univariate regression.

#' @param x A reg object
#' @param cuts A vector of significance values.
#' @param all_var A logical, whether to display all variables at specific p value cut-off, if \code{FALSE}, only display the first 10 variables
#' @param \dots additional arguments
#' @seealso \code{\link{dataframe.reg}},
#'   \code{\link{display}},  \code{\link{display.reg}}
#' @return Summary information of significant varibles.
#' @importFrom utils head
#' @export
#' @examples
#' reg_glm<-reg(data = diabetes, y = 5, factor = c(1, 3, 4), model = 'glm')
#' display(reg_glm)
#' display(reg_glm, all_var = FALSE)


display.reg <- function(x, cuts = c(0.001, 0.01, 0.05, 0.1, 1), all_var = TRUE,
    ...) {
    if (class(x) != "reg") {
        stop("x should be a `reg` object.", call. = FALSE)
    }
    cat("\nCall:\n", deparse(x$detail$call), "\n\n", sep = "")
    cat("Number of variables:", length(x$detail) - 1, "\n", sep = "\t")
    cat("Number of terms:", NROW(x$dataframe), "\n", sep = "\t")
    non_missing <- !is.na(x$dataframe$p.value)
    pvalue <- x$dataframe$p.value[non_missing]
    cat("Number of significant terms(alpha=0.05):", sum(pvalue < 0.05), "\n",
        sep = "\t")
    cat("\n")
    cat("Cumulative number of terms:\n")
    cat("\n")
    counts <- sapply(cuts, function(p) list(`p-value` = sum(pvalue <= p), var = x$dataframe[x$datafram$p.value <
        p, ]$term))


    counts <- t(counts)
    counts <- as.data.frame(counts)
    row.names(counts) <- paste("p < ", cuts, sep = "")
    names(counts)[1] <- "number of terms"
    print(counts[, 1, drop = FALSE])
    cat("\n\n")

    for (i in 1:NROW(counts)) {
        if (all_var) {
            cat(paste0(row.names(counts)[i], ":  ", paste(unlist(counts[i,
                2]), sep = ", ", collapse = ", ")), collapse = "\n\n\n")
        } else {
            cat(paste0(row.names(counts)[i], ":  ", paste(head(unlist(counts[i,
                2]), 10), sep = ", ", collapse = ", ")), collapse = "\n\n")
        }
    }

}
