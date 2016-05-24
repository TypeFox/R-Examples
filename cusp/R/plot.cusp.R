`plot.cusp` <-
function (x, what = c("all", "bifurcation", "residual", "densities"), 
    ...) 
{
    object = x
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    switch(match.arg(what, c("all", "bifurcation", "residual", 
        "densities")), bifurcation = plotCuspBifurcation(object, 
        ...), residual = plotCuspResidfitted(object, ...), 
        densities = plotCuspDensities(object, ...), all = {
            layout(matrix(c(1, 1, 1, 2, 1, 1, 1, 2, 2 + 1:4), 
                4))
            plotCuspBifurcation(object, ...)
            plotCuspResidfitted(object, ...)
            plotCuspDensities(object, ...)
        })
}

