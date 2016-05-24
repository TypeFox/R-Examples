`optimizeNetwork` <-
function (observations, predGrid, candidates, method, action,
    nDiff, model, criterion = "MUKV", plotOptim = TRUE, nGridCells,
    nTry, nr_iterations = 10000, formulaString, ...)
{
    if (length(method) == 0)
        stop(cat("No 'method' informed..."))
    if (length(method) > 0) {
        if (method != "spcov" & method != "ssa" & method != "manual")
            stop(cat("No 'method' required ..."))
        if (method == "ssa" & length(criterion) > 0) {
            if (criterion != "MUKV")
                stop(cat("Criterion ", criterion, " is not implemented."))
            if (criterion == "MUKV" & length(predGrid) == 0)
                stop(cat("No prediction locations to compute MUKV."))
        }
    }
    if (length(action) == 0)
        stop(cat("No 'action' defined ... choose 'add' or 'del'."))
    if (length(action) > 0) {
        if (action != "add" & action != "del")
            stop(cat("No relevant 'action' defined ... choose 'add' or 'del'."))
        if (action == "add") {
            if (length(candidates) == 0 | class(candidates) !=
                "SpatialPolygonsDataFrame")
                stop(cat("Candidate locations for additionnal measurements should be a shapefile."))
        }
    }
    if (length(nDiff) == 0 | nDiff <= 0)
        stop(cat("nDiff is not well defined", nDiff))
    if (method == "ssa") {
        if (missing(formulaString) || is.null(formulaString)) {
            observations = SpatialPointsDataFrame(observations,
                data = data.frame(dum = rep(1, dim(coordinates(observations))[1])))
            formulaString = dum ~ 1
        }
        return(ssaOptim(observations, predGrid, candidates, action,
            nDiff, model, nr_iterations, plotOptim, formulaString,
            ...))
    }
    if (method == "spcov") {
        if (action == "add") {
            return(spCovAdd(observations, candidates, nDiff,
                nGridCells, plotOptim, nTry))
        }
        if (action == "del") {
            return(spCovDel(observations, candidates, nDiff,
                plotOptim))
        }
    }
    if (method == "manual") {
        if (action == "add") {
            return(addManual(candidates, observations, nDiff))
        }
        if (action == "del") {
            return(delManual(candidates, observations, nDiff))
        }
    }
}
