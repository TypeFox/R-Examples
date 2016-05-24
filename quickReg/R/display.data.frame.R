#' Display summary information of varibles in a data.frame
#'
#' Dispaly count, frequency or mean, standard deviation and test of normality, etc.

#' @param x A data.frame
#' @param col Column indices of variables in the dataset to display, the default columns are all the variables
#' @param normtest  A character indicating test of normality, the default method is \code{\link{shapiro.test}} when sample size no more than 5000, otherwise \code{\link[nortest]{lillie.test}}{Kolmogorov-Smirnov} is used, see package \strong{nortest} for more methods.Use 'shapiro.test', 'lillie.test', 'ad.test', etc to specify methods.
#' @param useNA Whether to include NA values in the table, see \code{\link{table}} for more details
#' @param discrete_limit A numeric defining the minimal of unique value to display the variable as count and frequency
#' @param exclude_discrete Logical, whether to exclude discrete variables with more unique values specified by discrete_limit
#' @param \dots additional arguments
#' @import nortest
#' @export
#' @seealso \code{\link{display.reg}},  \code{\link{display}}
#' @examples
#' data(diabetes)
#' head(diabetes)
#' display(diabetes, 1:2)
#' display(diabetes, 1:11, normtest = "lillie.test")

display.data.frame <- function(x = NULL, col = NULL, normtest = NULL,  useNA = "ifany", discrete_limit = 10, exclude_discrete=TRUE, ...) {
    data <- x
    stopifnot(is.data.frame(data))
    if (is.null(col))
        col = seq_along(1:NCOL(data))
    if (is.null(normtest)) {
       normtest <- ifelse(NROW(data) <= 5000, "shapiro.test", "lillie.test")
    }
    result <- list()
    for (i in col) {
        split_line <- paste0(rep.int("=",80),collapse = "")
        term <- names(data)[i]
        if (length(unique(data[, i])) >=  discrete_limit && is.numeric(data[,
            i])) {
            summary <- summary(data[, i])
            describe <- psych::describe(data[, i])[2:13]
            row.names(describe)<-""

            normality <- do.call(normtest,list(data[, i]))
            digits <- getOption("digits")
            statistic<-format(signif(normality$statistic, max(1L, digits - 2L)))
            p.value<-format.pval(normality$p.value, digits = max(1L, digits - 3L))
            p.value<-paste("p-value", if (substr(p.value, 1L, 1L) == "<") p.value
                           else paste("=", p.value))
            result[[term]] <- list(split_line = split_line, summary = summary,
                describe = describe, normality = paste(normality$method,paste0("statistic = ", statistic), p.value,sep=", "))

        } else if (length(unique(data[, i])) >= discrete_limit&&exclude_discrete) {
          result[[term]] <- list(split_line = split_line, warning=paste0("`",term,"`"," is verbose. To display it using `exclude_discrete = FALSE` "))
          } else {
            count <- table(data[, i])
            propotion <- prop.table(count)

            if (useNA != "no") {
                if (any(is.na(data[, i]))) {
                  count_NA <- table(data[, i], useNA = useNA)
                  propotion_NA <- prop.table(count_NA)
                  result[[term]] <- list(split_line = split_line, table = rbind(count,
                    propotion), table_NA = rbind(count_NA, propotion_NA))
                } else {
                  result[[term]] <- list(split_line = split_line, table = rbind(count,
                    propotion))

                }
            }
        }
    }
    result
}
