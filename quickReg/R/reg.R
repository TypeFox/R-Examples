#' Apply univariate regression models
#'
#' Apply general linear model, generalized linear model, cox regression model,etc.

#' @param data A data.frame
#' @param x Column indices or names of the variables to be included in univariate analysis, the default columns are all the variables except dependent
#' @param y Column indice of dependent variable or name
#' @param factor Column indices of the variables or names to be treated as factor
#' @param model Univariate analysis method, see \code{\link{lm}}, \code{\link{glm}}, \code{\link[survival]{coxph}}
#' @param time Column indice of survival time or name, used in cox regression, see \code{\link[survival]{coxph}} for more details
#' @param \dots Further arguments passed to regression model
#' @return The return result is a list including two componets, the first part is a detailed anaysis result, the second part is a concentrated result in a  data.frame
#' @importFrom stats binomial confint glm lm
#' @export
#' @examples
#' reg_glm<-reg(data = diabetes, x = c(1:4, 6), y = 5, factor = c(1, 3, 4), model = 'glm')
#' ##  subset result like a list
#' reg_glm$detail
#' reg_glm$dataframe
#' reg_glm[2]
#' reg_glm$detail[2:4]
#' ##  other methods
#' reg(data = diabetes, x = c(1, 3:6), y = 10, factor = c(1, 3, 4), model = 'lm')
#' reg(data = diabetes, x = c( "sex","education","BMI"), y = "diabetes",
#' time ="age", factor = c("sex","smoking","education"), model = 'coxph')


reg <- function(data = NULL, x = NULL, y = NULL, factor = NULL, model = NULL,
    time = NULL, ...) {
  stopifnot(is.data.frame(data))
  if (is.null(y) || length(y) != 1)
    stop("One dependent varibale should be provided!", call. = FALSE)
  if (is.character(x)) x<-match(x,names(data))
  if (is.character(y)) y<-match(y,names(data))
  if (is.character(factor)) factor<-match(factor,names(data))
  if (is.character(time)) time<-match(time,names(data))
  if (is.null(x))
        x = setdiff(seq_along(1:NCOL(data)), c(y, time))
    if (y %in% x) {
      warning(paste0("Column indice ", y, " have overlap with x, please check it out."),
              call. = FALSE)
        x = setdiff(seq_along(1:NCOL(data)), c(y, time))
    }
    if (is.null(factor)) {
      data <- data
    } else  data[, factor] <- lapply(data[, factor], factor)

    if (!(model %in% c("lm", "glm", "coxph")))
        stop("model should be one of `lm`, `glm` or `coxph`.", call. = FALSE)
    arg <- list(...)
    if (any(c("data", "x", "y", "model") %in% names(arg)))
        stop("Please check the arguments: data, x, y, model!", call. = FALSE)

    result_detail <- result_dataframe <- list()

    split_line <- paste0(rep.int("=",80),collapse = "")
    for (i in x) {
        term <- names(data)[i]
        if (model == "lm") {
            fit <- lm(data[, y] ~ data[, i], ...)
            coef <- cbind(fit$coef, suppressMessages(confint(fit)))
            one <- cbind(term, summary(fit)$coefficients, coef)
            result_dataframe <- rbind(result_dataframe, one)
            result_detail[[term]] <- list(split_line = split_line, summary = summary(fit))
        } else if (model == "glm") {
            if ("family" %in% names(arg))
                fit <- glm(data[, y] ~ data[, i], arg) else {
                fit <- glm(data[, y] ~ data[, i], family = binomial(link = "logit"),
                  ...)
            }
            coef <- cbind(fit$coef, suppressMessages(confint(fit)))
            or = exp(coef)
            one <- cbind(term, summary(fit)$coefficients, or)
            result_dataframe <- rbind(result_dataframe, one)

            result_detail[[term]] <- list(split_line = split_line, summary = summary(fit),
                `OR(95%CI)` = or)
        } else if (model == "coxph") {
            fit <- survival::coxph(survival::Surv(data[, time], data[, y]) ~
                data[, i], ...)
            one <- cbind(term, summary(fit)$coefficients, exp(confint(fit)))
            result_dataframe <- rbind(result_dataframe, one)

            result_detail[[term]] <- list(split_line = split_line, summary = summary(fit))
        }

    }
    result_dataframe <- cbind(row.names = row.names(result_dataframe), result_dataframe)
    row.names(result_dataframe) <- NULL
    result_dataframe <- as.data.frame(result_dataframe, stringsAsFactors = FALSE)
    result_dataframe <- result_dataframe[result_dataframe$row.names != "(Intercept)",
        ]
    result_dataframe$row.names <- sub("data[, i]", "", result_dataframe$row.names,
        fixed = TRUE)
    result_dataframe$term <- paste0(result_dataframe$term, result_dataframe$row.names)
    result_dataframe$row.names <- NULL
    result_dataframe[, 2:8] <- lapply(result_dataframe[, 2:8], as.numeric)

    if (model == "lm") {
        names(result_dataframe) <- c("term", "estimate", "std.error", "statistic",
            "p.value", "coef", "coef.low", "coef.high")

    } else if (model == "glm") {
        names(result_dataframe) <- c("term", "estimate", "std.error", "statistic",
            "p.value", "OR", "OR.low", "OR.high")
    } else if (model == "coxph") {
        names(result_dataframe) <- c("term", "estimate", "HR", "std.error",
            "statistic", "p.value", "HR.low", "HR.high")
        result_dataframe <- result_dataframe[, c("term", "estimate", "std.error",
            "statistic", "p.value", "HR", "HR.low", "HR.high")]
    }

    result_detail[["call"]] <- match.call()
    result <- list(detail = result_detail, dataframe = result_dataframe)
    class(result) <- "reg"
    return(result)
}


