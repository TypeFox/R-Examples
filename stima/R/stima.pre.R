stima.pre <-
function (dataset, maxsplit, vfold = vfold, model="regtrunk", controlpre = NULL)
  { #default model is a regression trunk ; Classification trunk model is not available yet.  
REvec <- numeric(dim(dataset)[2] - 1)
    if (controlpre$sel == "backward") {
        controlpre$sel <- "none"
    }
    for (i in 2:dim(dataset)[2]) {
        obj <- stima(dataset, maxsplit, first = i, vfold = vfold,
            control = controlpre, printoutput = FALSE)
        if (maxsplit > dim(obj$goffull)[1]) {
            maxsplit <- dim(obj$goffull)[1]
        }
        if (vfold == 0) {
            REvec[i - 1] <- obj$goffull[maxsplit, 3]
        }
        else {
            REvec[i - 1] <-ifelse(maxsplit>2, min(obj$goffull[3:maxsplit, 5]), min(obj$goffull[2:maxsplit, 5]))
        }
    }
    best <- which(REvec == min(REvec)) + 1
    if (length(best) != 1) {
        cat("there are several best predictors for first split: ",
            names(dataset)[best], "\n")
        cat("we have selected the first one: ", names(dataset)[best][1],
            "\n")
        best <- best[1]
    }
    return(best)
}
