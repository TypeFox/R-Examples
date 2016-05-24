################################################################################
### Categorize functions and methods for a specific class
### (this is an internal utility function used in some of the package vignettes)
###
### Copyright (C) 2014-2016  Sebastian Meyer <Sebastian.Meyer@ifspm.uzh.ch>
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################


functionTable <- function (class, functions = list(),
                           format = "\\texttt", format.nongenerics = "\\textit")
{
    ## categorization of known generic functions
    KNOWNGENERICS <- list(
        Display = c("print", "summary", "xtable",
                    "plot", "animate", "as.stepfun",
                    "intensityplot"),
        Subset = c("[", "head", "tail", "subset"),
        Extract = c("nobs", "marks",
                    "coef", "fixef", "ranef", "vcov", "confint", "coeflist",
                    "logLik", "AIC", "extractAIC", "profile", "residuals",
                    "terms", "formula", "R0"),
        Modify = c("update", "untie", "add1", "drop1"),
        Convert = c("as.epidata"),
        Other = c("predict", "simulate", "pit", "scores", "calibrationTest")
    )
    
    if (is.null(names(functions)))  # put all functions in category "Other"
        functions <- list(Other = unlist(functions, use.names=FALSE))
    
    ## union known generics with specified functions
    categoryNames <- union(names(KNOWNGENERICS), names(functions))
    knowngenerics <- mapply(
        FUN = union, setNames(KNOWNGENERICS[categoryNames], categoryNames),
        functions[categoryNames], SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    ## get registered methods and associated generics
    allmethods <- methods(class = class)
    allgenerics <- attr(allmethods, "info")$generic
    genericsList <- lapply(X = knowngenerics, FUN = intersect, allgenerics)
    genericsList$Other <- c(genericsList$Other,
                            setdiff(allgenerics,
                                    unlist(genericsList, use.names=FALSE)))
    
    ## all extra 'functions' are not generic or without a method for 'class'
    nongenericsList <- lapply(X = functions, FUN = function (fnames) {
        res <- setdiff(fnames, allgenerics)
        ## note: we do not check if these functions actually exist()
        if (length(res)) paste0(format.nongenerics, "{", res, "}") else res
    })
    
    ## merge generics and non-generics
    functionList <- mapply(FUN = c, genericsList,
                           nongenericsList[names(genericsList)],
                           SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    ## transform list into a matrix by filling with empty cells
    categoryLengths <- lengths(functionList, use.names = FALSE)
    nrows <- max(categoryLengths)
    functionTable <- vapply(X = functionList[categoryLengths > 0L],
                            FUN = function (x)
                                c(paste0(format, "{", x, "}"),
                                  rep.int(NA_character_, nrows-length(x))),
                            FUN.VALUE = character(nrows),
                            USE.NAMES = TRUE)

    ## done
    functionTable #xtable::xtable(functionTable, ...)
}
