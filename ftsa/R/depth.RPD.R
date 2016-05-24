`depth.RPD` <- function (data, nproj = 50, deriv = c(0, 1), trim = 0.25, 
    dfunc2 = depth.mode, x = NULL, spline = TRUE, graph = FALSE, 
    ...) 
{
    functions = t(data$y)
    if (is.data.frame(functions)) 
        functions = as.matrix(functions)
    nr = nrow(functions)
    nc = ncol(functions)
    modulo = function(z) {
        sqrt(sum(z^2))
    }
    if (is.null(nr) || is.null(nc)) 
        stop("Input must be a matrix")
    if (is.null(x)) 
        x = 1:nc
    newfunc = array(NA, dim = c(nr, nc, length(deriv)))
    for (k in 1:length(deriv)) {
        if (deriv[k] == 0) {
            newfunc[, , k] = functions
        }
        else {
            for (i in 1:nr) {
                if (spline) {
                  fn = splinefun(x, functions[i, ])
                  newfunc[i, , k] = fn(x, deriv = deriv[k])
                }
                else {
                  newfunc[i, , k] = c(rep(0, deriv[k]), diff(functions[i, 
                    ], differences = deriv[k]))
                }
            }
        }
    }
    prof = rep(0, nr)
    vproject = matrix(0, nrow = nr, ncol = length(deriv))
    z = rnorm(nc * nproj)
    z = matrix(z, nrow = nproj, ncol = nc)
    modu = apply(z, 1, modulo)
    z = z/modu
    for (j in 1:nproj) {
        for (k in 1:length(deriv)) {
            matrx = newfunc[, , k]
            vproject[, k] = matrx %*% z[j, ]
        }
        resul = dfunc2(fts(1:dim(vproject)[2],t(vproject)), ...)
        prof = prof + resul$prof
    }
    prof = prof/nproj
    k = which.max(prof)
    med = functions[k, ]
    lista = which(prof >= quantile(prof, probs = trim))
    mtrim = apply(functions[lista, ], 2, mean)
    if (graph) {
        dev.new()
        cgray = 1 - (prof - min(prof))/(max(prof) - min(prof))
        if (nc == 2) {
            plot(range(functions[, 1]), range(functions[, 2]), 
                type = "n", xlab = colnames(functions)[1], ylab = colnames(functions)[2])
            points(functions[, 1], functions[, 2], col = gray(cgray))
            points(rbind(mtrim), pch = 19, col = gray(2 * trim), 
                cex = 2)
            points(rbind(med), col = 3, pch = 20, cex = 2)
        }
        else {
            plot(range(x), range(functions), type = "n", xlab = "t", 
                ylab = "X(t)", main = "RPD Depth")
            for (i in 1:nr) {
                lines(x, functions[i, ], col = gray(cgray[i]))
            }
            lines(x, mtrim, lwd = 2, col = "yellow")
            lines(x, med, col = "red", lwd = 2)
            legend("topleft", legend = c(paste("Trim", trim * 
                100, "%", sep = ""), "Median"), lwd = 2, col = c("yellow", 
                "red"))
        }
    }
    return(list(median = med, lmed = k, mtrim = mtrim, ltrim = lista, 
        prof = prof, proj = z))
}
