LKrig.make.grid.info<- function(grid.info, NC, x, V,distance.type, verbose=FALSE){
   if (is.null(grid.info)) {
        # if center grid information is missing create grid based on the locations
        # find range of scaled locations
        if (is.null(x)) {
            warning("Default spatial domain is  [0,1]X[0,1]")
            x<- cbind( c( 0,1), c(0,1))
        }
        if (is.null(NC)) {
            stop("need to specify NC for grid size")
        }
        range.x <- apply(as.matrix(x) %*% t(solve(V)), 2, "range")
        if (verbose) {
            cat("ranges of transformed variables", range.x, fill = TRUE)
        }
        grid.info <- list(xmin = range.x[1, 1], xmax = range.x[2, 
            1], ymin = range.x[1, 2], ymax = range.x[2, 2])
        d1 <- grid.info$xmax - grid.info$xmin
        d2 <- grid.info$ymax - grid.info$ymin
        grid.info$delta <- ifelse(distance.type == "cylinder", 
            d1/(NC - 1), max(c(d1, d2))/(NC - 1))
    }
    
    # check that cylinder geometry has the right delta with x i.e.divides range evenly.
    if (distance.type == "cylinder") {
        d1 <- grid.info$xmax - grid.info$xmin
        NC.test <- d1/grid.info$delta
        if (round(NC.test - round(NC.test)) > 1e-08) {
            stop("delta spacing in x dimension must be even for\ncylinder geometry")
        }
    }
     
# check distance choices
    if (is.na(match(distance.type, c("Euclidean", "cylinder")))) {
        stop("distance type is not supported (or is misspelled!).")
    }
#further checks for cylinder case
    if (distance.type == "cylinder") {
        if (V[1, 1] != 1) {
            stop("can not scale the angular coordinate (x[,1]),\nassumed to be in degrees")
        }
        if ((V[2, 1] != 0) | (V[1, 2] != 0)) {
            stop("can not have off diagonal elements in V")
        }
    }
   return(grid.info)
 }
