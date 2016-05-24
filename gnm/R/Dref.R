#  Copyright (C) 2005-2007, 2012 Heather Turner
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

Dref <- function(..., delta = ~ 1){
    preds <- match.call(expand.dots = FALSE)[["..."]]
    n <- length(preds)
    preds <- c(delta = rep(list(delta), n), preds)
    common <- c(1:n, rep(n + 1, n))
    extra <- setdiff(names(match.call()[-1]), c("", "delta"))
    if (length(extra))
        stop(paste(c("invalid argument passed to Dref:",
                     extra), collapse = " "))
    nf <- match(c("delta"), names(match.call()[-1]), 0)
    if ("formula" %in% names(match.call()[-1]))
        stop("formula argument of old plug-in has been renamed ",
             "\"delta\" in this function.")
    match <- c(rep(nf, n), 1:n)
    names(preds) <- c(rep("delta", n), rep("", n))

    list(predictors = preds,
         common = common,
         match = match,
         term = function(predLabels, ...){
             delta <- predLabels[1:n]
             gamma <- predLabels[-c(1:n)]
             paste("(((exp(", delta, "))/(",
                   paste("exp(", delta, ")", collapse = " + "),
                   "))*", gamma, ")", sep = "", collapse = " + ")},
         start = function(theta) {
             ifelse(attr(theta, "assign") == n + 1, 0.5,
                    runif(length(theta)) - 0.5)
         },
         call = as.expression(match.call()))
}
class(Dref) <- "nonlin"
