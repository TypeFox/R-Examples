## assumes convex objectives
`profZoom` <-
function (prof, ...) 
UseMethod("profZoom")
`profZoom.profileModel` <-
function (prof, max.zoom = 100, endpoint.tolerance = 0.001, verbose = FALSE, 
    ...) 
{
    if (is.null(prof$quantile)) 
        stop("An object with non-NULL quantile has to be supplied.")
    zero.bound <- prof$zero.bound
    intersects <- prof$intersects
    agreement <- prof$agreement
    fitted <- prof$fit
    which <- prof$profiled.parameters
    profRes <- prof$profiles
    isNA <- prof$isNA
    p <- length(profRes)
    if (scale <- !is.null(fitted$X.max.scaleFit)) {
        Xmax <- fitted$X.max.scaleFit
        for (i in 1:p) {
            if (isNA[i]) 
                next
            profRes[[i]][, 1] <- profRes[[i]][, 1] * Xmax[which[i]]
        }
    }
    betaNames <- names(profRes)
    quantile <- prof$quantile
    objective <- prof$profiled.objective
    grid.right <- matrix(Inf, p, 2)
    grid.left <- -grid.right
    for (i in 1:p) {
        if (isNA[i]) {
            grid.left[i, ] <- grid.right[i, ] <- NA
            next
        }
        profRes.i <- profRes[[i]]
        pos <- which(profRes.i[, 2] < quantile)
        if (intersects[i, 1]) 
            grid.left[i, ] <- profRes.i[pos[1] - c(1, 0), 1]
        if (intersects[i, 2]) 
            grid.right[i, ] <- profRes.i[pos[length(pos)] + c(0, 
                1), 1]
    }
    if (!max.zoom) {
        result <- cbind(rowSums(grid.left)/2, rowSums(grid.right)/2)
        dimnames(result) <- list(betaNames, c("Lower", "Upper"))
        return(result)
    }
    for (i in 1:p) {
        if (isNA[i]) {
            grid.left[i, ] <- grid.right[i, ] <- NA
            next
        }
        which.i <- which[i]
        bb <- coef(fitted)[which.i]
        if (verbose) 
            cat("Zooming for parameter", betaNames[i],"...\n")
        grid.left.i <- grid.left[i, ]
        grid.right.i <- grid.right[i, ]
        if (intersects[i, 1]) {
            zoom.step <- 1
            test <- TRUE
            while (zoom.step <= max.zoom & test) {
                mean.left <- mean(grid.left.i)
                profRes.left <- profiling(fitted, grid.bounds = c(mean.left, 
                  bb), gridsize = 5, verbose = FALSE, 
                  objective = objective, agreement = agreement, 
                  profTraces = FALSE, which = which.i, zero.bound = zero.bound)[[1]][1,]
                isLess.left <- profRes.left[2] < quantile
                if (isLess.left) 
                  grid.left.i <- c(grid.left.i[1], profRes.left[1])
                else grid.left.i <- c(profRes.left[1], grid.left.i[2])
                test <- abs(profRes.left[2] - quantile) > endpoint.tolerance
                zoom.step <- zoom.step + 1
            }
            #cat(" Left:",zoom.step)
            #cat(zoom.step,"left",i,bb,mean.left,"\n")
            #print(profRes.left[2],16)
            grid.left[i, ] <- grid.left.i
        }
        if (intersects[i, 2]) {
            zoom.step <- 1
            test <- TRUE
            while (zoom.step <= max.zoom & test) {
                mean.right <- mean(grid.right.i)
                profRes.right <- profiling(fitted, grid.bounds = c(bb, 
                  mean.right), gridsize = 5, verbose = FALSE, 
                  objective = objective, agreement = agreement, 
                  profTraces = FALSE, which = which.i, zero.bound = zero.bound)[[1]][5,]
                isLess.right <- profRes.right[2] < quantile
                if (isLess.right) 
                  grid.right.i <- c(profRes.right[1], grid.right.i[2])
                else grid.right.i <- c(grid.right.i[1], profRes.right[1])
                test <- abs(profRes.right[2] - quantile) > endpoint.tolerance
                zoom.step <- zoom.step + 1
            }
            #cat(" Right:",zoom.step,"\n")
            #cat(zoom.step,"right",i,bb,mean.right,"\n")
            #print(profRes.right[2],16)
            grid.right[i, ] <- grid.right.i
        }
    }
    result <- cbind(rowSums(grid.left)/2, rowSums(grid.right)/2)/(if (scale) 
        Xmax[which]
    else 1)
    dimnames(result) <- list(betaNames, c("Lower", "Upper"))
    result
}
