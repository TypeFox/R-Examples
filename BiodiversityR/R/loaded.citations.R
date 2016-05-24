`loaded.citations` <-
function() {
    loaded.packages <- (.packages())
    cat("Loaded packages: \n")
    cat(loaded.packages, "\n")
    cat("\n Citations: \n")
    standard.packages <- c("datasets", "grDevices", "graphics", "grid", "methods", "splines", 
        "stats", "stats4", "tcltk", "tools", "utils")
    non.standard <- as.logical(rep("T", l=length(loaded.packages)) )
    for (i in 1:length(non.standard)) {
        if (any (loaded.packages[i] == standard.packages)) {non.standard[i] <- F}
    }
    loaded.packages <- loaded.packages[non.standard]
    for (i in 1:length(loaded.packages)) {print(utils::citation(loaded.packages[i]))}
}

