`depth.RPD` <- function (data, nproj = 50, deriv = c(0, 1), trim = 0.25, 
    dfunc2 = depth.RP, x = NULL, spline = TRUE, ...) 
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
        resul = dfunc2(fts(1:dim(vproject)[2], t(vproject)), trim = trim, ...)
        prof = prof + resul$prof
    }
    prof = prof/nproj
    k = which.max(prof)
    med = functions[k, ]
    lista = which(prof >= quantile(prof, probs = trim))
    mtrim = apply(functions[lista, ], 2, mean)
    return(list(median = med, lmed = k, mtrim = mtrim, ltrim = lista, 
        prof = prof))
}

