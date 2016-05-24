# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# fit rank-deficient generalized linear model
rdglm.fit <- function(x, y, weights = rep(1, nobs), start = NULL,
                      etastart = NULL, mustart = NULL, offset = rep(0, nobs),
                      family = gaussian(), control = list(), 
                      method = "firthglm.fit", intercept = TRUE)
{
    # method
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")

    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    nobs <- NROW(y)
    nvars <- ncol(x)

    # handle degenerate (no-observations) case
    if (nobs == 0L) {
        coef <- numeric(nvars)
        names(coef) <- xnames

        if (family$family %in% c("poisson", "binomial")) {
            dispersion <- 1
        } else {
            dispersion <- 0
        }

        return(list(coefficients = coef, rank = 0L, qr = qr(x),
                    df.residual = 0L, dispersion = dispersion,
                    prior.weights = weights))
    }

    fit <- eval(call(if (is.function(method)) "method" else method,
        x = x, y = y, weights = weights, start = start, etastart = etastart,
        mustart = mustart, offset = offset, family = family, control = control,
        intercept = intercept))

    coef <- fit$coefficients
    coef[is.na(coef)] <- 0

    df.residual <- fit$df.residual
    if (family$family %in% c("poisson", "binomial")) {
        dispersion <- 1
    } else if (df.residual > 0) {
        dispersion <- (sum((fit$weights * fit$residuals^2)[fit$weights > 0])
                       / df.residual)
    } else {
        dispersion <- 0
    }

    list(coefficients = coef, rank = fit$rank, qr = fit$qr,
         df.residual = df.residual, dispersion = dispersion,
         prior.weights = fit$prior.weights)
}
