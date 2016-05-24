krig <-
function(values, coords, grid, gv, m=NA, cv=FALSE, clamp=FALSE, verbose=TRUE) {

    if (!is.vector(values)) stop("Values must be a vector.")

    if (!class(grid) == "data.frame") {
        ## Attempt coercion to data.frame
        grid <- as.data.frame(grid)
    }
    
    n <- length(values) # number of samples
    ni <- nrow(grid)   # number of locations to interpolate

    # more error checking
    if (n != nrow(coords)) {
        stop('Number of values different from number of locations!')
    }

    mtest <- mtest.gv(gv)
    if (!mtest) stop("Object has no model built!")

    sill <- gv$model$sill

    if(verbose)
        if (is.na(m)) 
            print("Ordinary Kriging") else print("Simple Kriging")

    # Calculates simple Euclidean distances
    C <- as.matrix(dist(coords))
    Cc <- sill - predict(gv, C)
    if (is.na(m)) Cc <- rbind(cbind(Cc, rep(1, n)), c(rep(1, n), 0))
    #iC <- solve(Cc) # Inverse of the  matrix
    iC <- mpinv(Cc)

    Z <- sd <- rep(NA, ni)
    for (i in 1:ni) {
        pnt <- grid[i,]
        i.dists <- ((coords[,1]-pnt[,1])**2+(coords[,2]-pnt[,2])**2)**0.5
        D <- sill - predict(gv, i.dists)
        if (is.na(m)) D <- c(D, 1)
        w <- iC %*% D
        if (is.na(m))
            Z[i] <- sum(values*w[1:n, 1]) else
            Z[i] <- m + sum((values-m)*w[,1])
        sd[i] <- (sill - sum(w[,1]*D))**0.5
    }

    if (clamp) {
        Z[Z > 1] = 1
        Z[Z < 0] = 0
    }

    if (cv) {
        # Perform cross validation and mean squared error
        cv <- matrix(NA, length(values), 2)
        colnames(cv) <- c('value', 'predicted')
        for (i in 1:length(values)) {
            cv.k <- krig(values[-i], coords[-i,], coords[i,], 
                         gv=gv, m=m, cv=FALSE, verbose=FALSE)
            cv[i,] <- c(values[i], cv.k[1,1])
        }
        MSE <- sum(apply(cv, 1, diff)**2)/length(values)

        return(list(intpl=data.frame(Z, sd), cross.valid=cv, MSE=MSE))
    
    } else {
        # Return only interpolation
        return(data.frame(Z, sd))
    }
}


