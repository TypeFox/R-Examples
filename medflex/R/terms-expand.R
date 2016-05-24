#' Expanded dataset
#'
#' @description Expanded dataset including either ratio-of-mediator probability weights or imputed nested counterfactual outcomes.
#' @return A data frame, resulting from applying \code{\link{neWeight}} or \code{\link{neImpute}} on an original dataset \code{data}.
#' This data frame has \code{nRep * length(data)} rows, containing all original variables (except the original exposure variable)
#' and two variables reflecting observed and hypothetical values of the exposure for each observation unit.
#'
#' These auxiliary variables (\emph{x} and \emph{x*}) are named after the exposure variable and carry integers as suffixes.
#' Suffixes \code{0} and \code{1} are used for variables whose corresponding parameters in the final natural effect model index natural direct and indirect effects, respectively.
#'
#' This object also stores some additional attributes, which are used as input for \code{\link{neModel}}, such as
#'   \item{\code{model}}{the fitted working model object}
#'   \item{\code{data}}{original dataset}
#'   \item{\code{call}}{the matched call}
#   \item{\code{terms}}{the \code{\link{neTerms}} object used}
#'   \item{\code{terms}}{the \code{neTerms} (internal class) object used}
#'   \item{\code{weights}}{ratio-of-mediator probability weights (only stored if object inherits from class \code{weightData})}
#' @note If the weighting-based approach (\code{\link{neWeight}}) is applied, the original outcome values are copied for the nested counterfactual outcomes
#' and the object stores an additional attribute, \code{"weights"}, containing a vector with ratio-of-mediator probability weights.
#'
#' If the imputation-based approach (\code{\link{neImpute}}) is applied, the nested counterfactual outcomes are imputed by predictions from the imputation model.
#'
#' In the former case, this object inherits from classes \code{c("data.frame", "expData", "impData")}, whereas in the latter case it inherits from classes \code{c("data.frame", "expData", "weightData")}.
#' @seealso \code{\link{neImpute}}, \code{\link{neWeight}}
#' @name expData
NULL

expandData <- function (x, data, nMed, ...) 
{
    args <- eval(substitute(alist(x, data, ...)))
    if (inherits(data, "environment")) data <- as.data.frame(as.list(data))
    args[[1]] <- if (!is.null(args$vartype) && grepl("factor", 
        attr(eval(args$vartype), "xasis"))) {
        quote(as.factor(data[, x]))
    }
    else {
        quote(data[, x])
    }
    Xexp <- do.call("expandX", args)
    x <- attr(Xexp, "x")
    nRep <- ifelse(is.factor(x), nlevels(x), args$nRep)
    joint <- TRUE
    nExp <- 1
    if (nExp > 1) {
        ids <- matrix(seq.int(nrow(Xexp)), nrow = nrow(Xexp)/nRep, 
            byrow = TRUE)
        ids <- matrix(rep(ids, times = nMed), nrow = nrow(Xexp)/nRep)
        ids <- as.vector(t(ids))
        tmp <- rep(seq.int(nMed), each = nRep)
        tmp <- rep(tmp, times = length(x))
        Xexp <- cbind(tmp, Xexp[ids, ])
        Xexp <- cbind(Xexp, mapply(function(x) ifelse(Xexp[, 
            "tmp"] <= x - 1, Xexp[, "aux0"], Xexp[, "aux1"]), 
            2:nMed))
        if (is.factor(x)) 
            Xexp[, ncol(Xexp) - seq(0, nExp - 2)] <- apply(Xexp[, 
                ncol(Xexp) - seq(0, nExp - 2)], 2, function(y) factor(y, 
                labels = levels(Xexp[, "aux0"])))
        colnames(Xexp) <- c(colnames(Xexp)[1:3], paste0("aux", 
            2:nMed))
        Xexp <- Xexp[, -1]
    }
    data$id <- rownames(data)
    Xexp$id <- rep(rownames(data), each = nRep * nExp)
    expData <- suppressWarnings(merge(data, Xexp, by = "id", 
        sort = FALSE))
    return(expData)
}

expandX <- function (x, data, ...) 
{
    UseMethod("expandX") 
}

expandX.factor <- function (x, data, ...) 
{
    args <- as.list(match.call())[-1L]
    aux0 <- rep(x, each = nlevels(x))
    tmp1 <- as.numeric(x)
    aux1 <- factor(as.vector(mapply(function(y) (y + seq.int(nlevels(x)) - 
        2) %% nlevels(x) + 1, tmp1)), labels = levels(x))
    res <- data.frame(aux0, aux1)
    attr(res, "x") <- x
    return(res)
}

expandX.numeric <- function (x, data, nRep, xSampling, xFit, percLim, vartype, ...) 
{
    aux0 <- rep(x, each = nRep)
    if (missing(xFit)) {
        tmp <- if (length(vartype$C) == 0) 
            "1"
        else paste(vartype$C, collapse = "+")
        xFit <- as.formula(paste(vartype$X, tmp, sep = "~"))
        Xfit <- glm(xFit, data = data)
    }
    else {
        Xfit <- update(xFit, data = data)
    }
    Xmean <- predict.glm(Xfit, type = "response")
    Xsd <- with(summary.glm(Xfit), sqrt(deviance/df.residual))
    switch(xSampling, quantiles = {
        p <- if (is.null(percLim)) seq(0, 1, length.out = nRep + 
            2)[-c(1, nRep + 2)] else seq(percLim[1], percLim[2], 
            length.out = nRep)
        aux1 <- as.vector(outer(p, Xmean, function(p, mean) qnorm(p, 
            mean, sd = Xsd)))
    }, random = {
        aux1 <- as.vector(outer(rep.int(1, nRep), Xmean, function(n, 
            mean, sd) rnorm(n, mean, sd = Xsd)))
    })
    res <- as.data.frame(cbind(aux0, aux1))
    attr(res, "x") <- x
    return(res)
}

neTerms <- function (x, ...)
{
    UseMethod("neTerms")
}

neTerms.Mformula <- function (x, Y, ...) 
{
    nMed <- 1
    terms <- terms.formula(as.formula(x))
    predvars <- all.vars(terms[[3]])
    attr(terms, "vartype") <- list(Y = Y, X = predvars[1], M = c(if (nMed > 
        1) predvars[2:nMed], all.vars(terms[[2]])), C = predvars[-(1:nMed)])
    attr(attr(terms, "vartype"), "xasis") <- dimnames(attr(terms, 
        "factors"))[[1]][rowSums(attr(terms, "factors")) > 0][1]
    class(terms) <- c("neTerms.object", class(terms))
    return(terms)
}

neTerms.Xformula <- function (x, ...) 
{
    terms <- terms.formula(as.formula(x))
    predvars <- all.vars(terms[[3]])
    attr(terms, "vartype") <- list(X = all.vars(terms[[2]]), 
        C = predvars)
    class(terms) <- c("neTerms.object", class(terms))
    return(terms)
}

neTerms.Yformula <- function (x, nMed, ...) 
{
    terms <- terms.formula(as.formula(x))
    predvars <- all.vars(terms[[3]])
    attr(terms, "vartype") <- list(Y = all.vars(terms[[2]]), 
        X = predvars[1], M = predvars[2:(nMed + 1)], C = predvars[-(1:(nMed + 
            1))])
    attr(attr(terms, "vartype"), "xasis") <- dimnames(attr(terms, 
        "factors"))[[1]][rowSums(attr(terms, "factors")) > 0][1]
    attr(attr(terms, "vartype"), "masis") <- dimnames(attr(terms, 
        "factors"))[[1]][rowSums(attr(terms, "factors")) > 0][2:(nMed + 1)]
    class(terms) <- c("neTerms.object", class(terms))
    return(terms)
}
