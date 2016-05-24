#  Copyright (C) 2006, 2009, 2013, 2014 Heather Turner
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

expandCategorical <- function(data, catvar, sep = ".", countvar = "count",
                              idvar = "id", as.ordered = FALSE, group = TRUE) {

    cat <- interaction(data[catvar], sep = sep)
    ncat <- nlevels(cat)

    covar <- data[, -match(catvar, names(data)), drop = FALSE]
    catvar <- paste(catvar, collapse = sep)

    if (group == TRUE) {
        if (length(covar)) {
            ord <- do.call("order", covar)
            vars <- covar[ord, , drop = FALSE]
            dupvars <- duplicated(vars)
            d <- diff(c(which(!dupvars), length(dupvars) + 1))
            n <- sum(!dupvars)
            id <- factor(rep(seq(n), d))
            counts <- as.data.frame(table(list(cat = cat[ord], id = id)))
            newData <- vars[which(!dupvars)[counts$id], , drop = FALSE]
            rownames(newData) <- NULL
            newData[c(catvar, idvar, countvar)] <- counts
        } else {
            newData <- data.frame(table(cat))
            colnames(newData) <- c(catvar, countvar)
            newData[[idvar]] <- factor(1)
        }
    }
    else {
        n <- nrow(covar)
        id <- gl(n, ncat)
        newData <- covar[id, , drop = FALSE]
        newData[[catvar]] <-
            gl(ncat, 1, n * ncat, labels = levels(cat), ordered = as.ordered)
        newData[[countvar]] <- as.vector(t(class.ind(cat)))
        newData[[idvar]] <- id
    }
    newData
}
