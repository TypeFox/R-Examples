handlesPointAtEdge <-
function (point, fullerhull, fixedlist)
{
    if (length(fullerhull$b) > 2L) {
        fullpoint <- rep(NA, ncol(fullerhull$a))
        names(fullpoint) <- colnames(fullerhull$a)
        fullpoint[names(point)] <- point
        fullpoint[names(fixedlist)] <- unlist(fixedlist)
        mg <- fullerhull$a %*% fullpoint - fullerhull$b
        focalConstraints <- which(sapply(mg, function(v) {
            v == 0
        }))
        if (length(focalConstraints) > 0) {
            which.maskingpts <- which(unlist(lapply(fullerhull$affectedConstraints,
                function(v) {
                  length(focalConstraints %w/o% v) == 0
                })))
        }
        else which.maskingpts <- NULL
        masksize <- length(which.maskingpts)
        if (masksize == 0) {
            return(list(edgelevel = 0L, edge = NULL))
        }
        else if (masksize == 1) {
            edge <- fullerhull$vertices[which.maskingpts, , drop = F]
            return(list(edgelevel = 1L, edge = edge))
        }
        else if (masksize == 2L) {
            edge <- fullerhull$vertices[which.maskingpts, ]
            resu <- optimize(function(r) {
                tofKpredict.nohull(r * edge[1L, ] + (1 - r) *
                  edge[2L, ], fixedlist = fixedlist)
            }, lower = 0, upper = 1, maximum = T)
            rr <- resu$maximum
            resu$par <- rr * edge[1, ] + (1 - rr) * edge[2, ]
            names(resu)[which(names(resu) == "objective")] <- "value"
            resu <- c(resu, convergence = 0)
            return(list(edgelevel = 2, edge = edge, resu = resu))
        }
        else {
            edge <- t(t(fullerhull$vertices[which.maskingpts, ]))
            ui <- rbind(-1, diag(masksize - 1))
            ci <- c(-1, rep(0, masksize - 1))
            objfn <- function(rv) {
                rv <- c(rv, 1 - sum(rv))
                tofKpredict.nohull((rv %*% edge)[1L, ], fixedlist = fixedlist)
            }
            objfn.grad <- function(rv) {
                grad(func = objfn, x = rv)
            }
            resu <- constrOptim(rep(1/masksize, masksize - 1),
                objfn, objfn.grad, ui = ui, ci = ci, mu = 1e-08,
                method = "BFGS", control = list(fnscale = -1,
                  trace = FALSE, maxit = 10000))
            resu$par <- (c(resu$par, 1 - sum(resu$par)) %*% edge)[1,
                ]
            return(list(edgelevel = masksize, edge = edge, resu = resu))
        }
    }
    else return(list(edgelevel = 555, edge = NULL))
}
