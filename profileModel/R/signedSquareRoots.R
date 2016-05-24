`signedSquareRoots` <-
function (prof) 
UseMethod("signedSquareRoots")
`signedSquareRoots.profileModel` <-
function (prof) 
{
    if (!prof$agreement) 
        stop("The objective and the fitting procedure ", fitted$call[[1]], 
            " do not agree. Signed square roots cannot be calculated.")
    which <- prof$profiled.parameters
    isNA <- prof$isNA
    intersects <- prof$intersects
    fit <- prof$fit
    beta <- fit$coefficients[which]
    profRes <- prof$profiles
    p <- length(profRes)
    for (i in 1:p) {
        if (isNA[i]) 
            next
        profRes.i <- profRes[[i]]
        sgn <- sign(profRes.i[, 1] - beta[i])
        if (!is.null(intersects)) 
            if (sum(intersects[i, ]) == 1) 
                sgn <- sum(c(-1, 1) * intersects[i, ])
        sgn.sqrt <- sgn * sqrt(profRes.i[, 2])
        profRes[[i]][, 2] <- sgn.sqrt
        colnames(profRes[[i]])[2] <- "Signed sqrt of the objective"
    }
    prof$profiles <- profRes
    prof
}
