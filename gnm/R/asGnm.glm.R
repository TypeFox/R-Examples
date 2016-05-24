#  Copyright (C) 2006, 2008, 2010 Heather Turner
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

asGnm.glm <- function(object, ...) {
    glmExtra <- match(c("effects", "R", "qr", "null.deviance", "df.null",
                        "boundary", "control", "contrasts"), names(object))
    modelData <- model.frame(object)
    object[glmExtra] <- NULL
    object$call[[1]] <- as.name("gnm")
    constrain <- which(is.na(coef(object)))
    object$terms <- gnmTerms(object$formula, data = modelData)
    object <- c(list(eliminate = NULL, ofInterest = NULL,
                     na.action = na.action(modelData), constrain = constrain,
                     constrainTo = numeric(length(constrain))),
                object)
    names(object)[match("linear.predictors", names(object))] <- "predictors"
    if (is.null(object$offset))
        object$offset <- rep.int(0, length(coef(object)))
    object$tolerance <- object$iterStart <- object$iterMax <-
        "Not available - model fitted by glm()"
    class(object) <- c("gnm", "glm", "lm")
    object
}
