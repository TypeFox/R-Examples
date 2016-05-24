plot.vhica <-
function (x, sp1 = NULL, sp2 = NULL, ...) 
{
    op <- par(no.readonly = TRUE)
    species <- c(sp1, sp2)
    if (is.null(species)) {
        species <- unique(c(x$div[, ncol(x$div) - 0:1]))
    }
    if (length(species) < 2) 
        stop("VHICA analysis requires at least two species")
    ly <- matrix(0, ncol = length(species) - 1, nrow = length(species) - 
        1)
    ly[upper.tri(ly, diag = TRUE)] <- 1:(length(species) * (length(species) - 
        1)/2)
    layout(ly)
    for (sp1 in 1:(length(species) - 1)) {
        for (sp2 in (sp1 + 1):length(species)) {
            cross <- paste(species[sp1], species[sp2], sep = "X")
            if (!cross %in% names(x$reg)) {
                cross <- paste(species[sp2], species[sp1], sep = "X")
            }
            if (!cross %in% names(x$reg)) {
                stop("The cross between species ", species[sp1], 
                  " and ", species[sp2], " is not documented")
            }
            .plot.regression(x$reg[[cross]], ...)
        }
    }
    par(op)
}
