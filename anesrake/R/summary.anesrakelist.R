summary.anesrakelist <-
function(object, ...) {
    namer <- c("convergence", "raking.variables", "weight.summary")
    part1list <- list(paste(object$converge, "after", object$iterations, 
        "iterations"), object$varsused, summary(object$weightvec))
    names(part1list) <- namer
    part2list <- weightassess(object$targets, object$dataframe, 
        object$weightvec, object$prevec)
    out <- c(part1list, part2list)
    out
}

