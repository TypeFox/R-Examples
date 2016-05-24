#  Copyright (C) 2005, 2008, 2010, 2012 Heather Turner
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

termPredictors.gnm <- function(object, ...) {
    if (is.null(object$termPredictors)){
        modelData <- model.frame(object)
        modelTerms <- terms(object)
        if (!is.empty.model(modelTerms)) {
            modelTools <- gnmTools(modelTerms, modelData)
            theta <- parameters(object)
            varPredictors <- modelTools$varPredictors(theta)
            termPredictors <- modelTools$predictor(varPredictors, term = TRUE)
            rownames(termPredictors) <- rownames(modelData)
        }
        else termPredictors <- modelData[,0]
        if (!is.null(object$eliminate)) termPredictors <-
            cbind("(eliminate)" =
                  as.vector(attr(coef(object), "eliminated")[object$eliminate]),
                  termPredictors)
        termPredictors
    }
    else object$termPredictors
}
