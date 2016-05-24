"plt.cc" <-
function (res, d1 = 1, d2 = 2, int = 0.5, type = "b", ind.names=NULL, var.label=FALSE, Xnames=NULL, Ynames=NULL) 
{
    par(mfrow = c(1, 1), pty = "s")
    if (type == "v") 
        plt.var(res, d1, d2, int, var.label, Xnames, Ynames)
    if (type == "i") 
        plt.indiv(res, d1, d2, ind.names)
    if (type == "b") {
        def.par <- par(no.readonly = TRUE)
        layout(matrix(c(0, 0, 1, 2, 0, 0), ncol = 2, nrow = 3, 
        byrow = TRUE), widths = 1, heights = c(0.1, 1, 0.1))
        par(pty = "s", mar = c(4, 4.5, 0, 1))
        plt.var(res, d1, d2, int, var.label, Xnames, Ynames)
        plt.indiv(res, d1, d2, ind.names)
        par(def.par)
    }
}

