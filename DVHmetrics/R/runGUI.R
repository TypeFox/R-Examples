runGUI <-
function(...) {
    appDir <- system.file("DVHshiny", package="DVHmetrics")
    if(!nzchar(appDir)) {
        stop("Could not find Shiny directory. Try re-installing 'DVHmetrics'.", call.=FALSE)
    }

    shiny::runApp(appDir, ...)    
}
