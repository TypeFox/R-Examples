#  Modification of anova.glm from the stats package for R.
#
#  Copyright (C) 1995-2005 The R Core Team
#  Copyright (C) 2005, 2006, 2008, 2012 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

anova.gnm <- function (object, ..., dispersion = NULL, test = NULL)
{
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named))
        warning("the following arguments to 'anova.gnm' are invalid and dropped: ",
                paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.gnm <- unlist(lapply(dotargs, function(x) inherits(x, c("gnm", "glm"))))
    dotargs <- dotargs[is.gnm]
    if (length(dotargs) > 0)
        return(anova(structure(c(list(object), dotargs), class="glmlist"),
                     dispersion = dispersion, test = test))

    x <- model.matrix(object)
    varlist <- attr(terms(object), "term.labels")
    varseq <- attr(x, "assign")
    pars <- setdiff(unique(varseq), c(0, varseq[object$constrain]))

    nvars <- length(varlist)

    nonlinear <- match(TRUE, attr(terms(object), "type") != "Linear")

    if (is.na(nonlinear))
        nonlinear <- nvars + 1

    resdev <- resdf <- fit <- NULL
    origConstrain <- object$constrain
    origConstrainTo <- object$constrainTo
    if (nvars > 0) {
        for (i in pars) {
            if (i < nonlinear && is.null(object$eliminate)){
                fit <- glm.fit(x = x[, varseq < i, drop = FALSE],
                               y = c(object$y), offset = c(object$offset),
                               start = object$start,
                               weights = c(object$prior.weights),
                               family = object$family)
            }
            else {
                f <- update.formula(formula(object),
                                    paste(". ~ . -",
                                          paste(varlist[i:nvars],
                                                collapse = " - ")))
                f <- update.formula(formula(object), f)
                fit <- update(object, formula = f, verbose = FALSE)
            }
            resdev <- c(resdev, fit$deviance)
            resdf <- c(resdf, fit$df.residual)
        }
        resdf <- c(resdf, object$df.residual)
        resdev <- c(resdev, object$deviance)
        table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))),
                            resdf, resdev)
    }
    else
        table <- data.frame(NA, NA, object$df.residual, object$deviance)
    dimnames(table) <- list(c("NULL", labels(object)),
                            c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
    title <- paste("Analysis of Deviance Table", "\n\nModel: ",
                   object$family$family, ", link: ", object$family$link,
                   "\n\nResponse: ", as.character(formula(object)[[2]]),
                   "\n\nTerms added sequentially (first to last)\n\n",
                   sep = "")
    df.dispersion <- Inf
    if (is.null(dispersion)) {
        dispersion <- attr(vcov(object), "dispersion")
        df.dispersion <- if (dispersion == 1)
            Inf
        else object$df.residual
    }
    if (!is.null(test))
        table <- stat.anova(table = table, test = test, scale = dispersion,
            df.scale = df.dispersion, n = NROW(x))
    structure(table, heading = title, class = c("anova", "data.frame"))
}
