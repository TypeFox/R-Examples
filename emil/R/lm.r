#' Fit a linear model fitted with ordinary least squares
#' 
#' Based on \code{\link{lm}}.
#' 
#' @param x Descriptors.
#' @param y Response, numeric.
#' @param formula See \code{\link{lm}}.
#' @param ... Sent to \code{\link{lm}}.
#' @return Fitted linear model.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{predict_lm}},
#'   \code{\link{modeling_procedure}}
#' @export
fit_lm <- function(x, y, formula=y~., ...){
    df <- data.frame(y, x)
    rm(y,x)
    vars.missing <- setdiff(all.vars(formula), c(".", names(df)))
    if(!is_blank(vars.missing)){
        omitted <- length(vars.missing) - 20
        vars.missing <- paste(sprintf("`%s`", head(vars.missing, 20)), collapse=", ")
        if(omitted > 0) vars.missing <- paste(vars.missing, "+", omitted, "more")

        omitted <- length(names(df)) - 20
        vars.present <- paste(sprintf("`%s`", head(names(df), 20)), collapse=", ")
        if(omitted > 0) vars.present <- paste(vars.present, "+", omitted, "more")

        stop(sprintf("Variables %s not found in data frame. Variables available are %s.",
            vars.missing, vars.present))
    }
    lm(formula, data.frame(y, x), ...)
}

#' Prediction using linear model
#' 
#' @param object Fitted classifier produced by \code{\link{fit_lm}}.
#' @param x Dataset to be predicted upon.
#' @param ... Sent to \code{\link{predict.lm}}
#' @return A list with elements:
#' \itemize{
#'     \item{\code{prediction}: Vector of predicted response.}
#' }
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{emil}}, \code{\link{fit_lm}},
#'   \code{\link{modeling_procedure}}
#' @export
predict_lm <- function(object, x, ...){
    list(prediction = predict.lm(object, data.frame(y=NA, x), ...))
}

