runGUI <-
function(app=c("analyze", "hitprob", "range", "angular"), ...) {
    app <- match.arg(toupper(app),
                     choices=c("ANALYZE", "HITPROB", "RANGE", "ANGULAR"),
                     several.ok=FALSE)

    appDir <- if(app == "ANALYZE") {
        system.file("shinyAnalyzeGroups", package="shotGroups")
    } else if(app == "HITPROB") {
        system.file("shinyHitProbability", package="shotGroups")
    } else if(app == "RANGE") {
        system.file("shinyRangeStatistics", package="shotGroups")
    } else if(app == "ANGULAR") {
        system.file("shinyAngularSize", package="shotGroups")
    }

    if(!nzchar(appDir)) {
        stop("Could not find Shiny directory. Try re-installing 'shotGroups'.", call.=FALSE)
    }

    if(requireNamespace("shiny", quietly=TRUE)) {
        shiny::runApp(appDir, ...)
    } else {
        stop("Could not find package shiny - please install first")
    }
}
