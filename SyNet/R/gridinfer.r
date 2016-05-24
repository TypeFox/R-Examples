gridinfer <- function (file = NULL, dntable = NULL, sp_row = TRUE, reciprocity = TRUE,
                       criterion = "max", tolerance = sqrt(2), conditioned = TRUE, ...)
{
    coords <- NULL
    if (is.null(dntable))
        dntable <- read.table(file, ...)
    if(sp_row) m <- as.matrix(dntable) else m <- as.matrix(t(dntable))
    if(any(m < 0)) stop("Sorry, but distributional matrix includes misleading negative entries\n")
    coords <- floor(t(m[1:2,,drop=FALSE]))
    if(any(duplicated(coords))) stop(paste("Each cell must be identified through an exclusive pair of coordinates\n",
                                     "Please, check your coordinates for duplciates"))
    m <- m[-(1:2),,drop = FALSE] & m[-(1:2),,drop =FALSE]
    if(!all(apply(m, 1, any))) stop("Error: you have included species without presences")
    dcells <- as.matrix(dist(coords))
    if(conditioned) sm <- m %*% t(m) else sm <- matrix(TRUE, nrow(m), nrow(m))
    if(reciprocity) oper <- "&" else oper <- "|"
    for(i in 1:(nrow(m)-1))
      for(j in (i+1):nrow(m)) {
         sp1 <- apply(dcells[m[i,], m[j,], drop = FALSE], 1, min)
         sp2 <- apply(dcells[m[i,], m[j,], drop = FALSE], 2, min)
         if(!(eval(call(oper, eval(call(criterion, sp1)) <= tolerance,
         eval(call(criterion, sp2)) <= tolerance)) & sm[i, j])) sm[i,j] <- sm[j, i] <- FALSE
    }
    out <- list(sm = ifelse(sm, 1, 0), Label = rownames(m), occupancy = apply(m, 1, which), coords = coords, kind = "grids")
    class(out) <- "gridinference"
    return(out)
}

