assimilation <-
function(plants, # marked point pattern (class ppp)
         pixsize=0.2, # resolution (pixel size)
         resource=1, # resource map
         influence=gnomon.inf, # influence function
         infpar=list(a=1, b=4, smark=1), # influence function parameter(s)
         asym=Inf, # asymmetry parameter
         efficiency=flat.eff, # efficiency function
         effpar=NULL, # efficiency function parameter(s)
         plot=TRUE, # plot influences map?
         afree=FALSE, # include free-growing assimilation in output?         
         centroid=FALSE # include assimilation centroids in output?
        )
{
    # Resource pixel image
    resource <- as.im(resource, as.mask(plants, eps=pixsize))
    resmatrix <- as.matrix(resource)
    # This is dangerous, spatstat warns about possible implementation changes:
    xcol <- resource$xcol
    yrow <- resource$yrow
    pixarea <- resource$xstep * resource$ystep
    # Extract info
    xplant <- coords(plants)$x
    yplant <- coords(plants)$y
    marks <- as.data.frame(marks(plants))
    # First pass, acccumulate partition function denominator
    infi <- function(y, x) influence(dx = x - xplant[i], dy = y - yplant[i],
        marks = marks[i, ], par = infpar)
    denom <- 0 * resmatrix
    for(i in 1:npoints(plants)) {
        phi <- outer(yrow, xcol, infi)
        if(is.infinite(asym)) denom <- pmax(denom, phi)
        else if(asym == 1) denom <- denom + phi
        else if(asym > 0) denom <- denom + phi^asym
        else denom <- denom + (phi > 0) # asym == 0
    }
    if(plot) plot(as.im(denom, resource), main=NULL)
    # Second pass, get assimilation index for each plant
    effi <- function(y, x) efficiency(dx = x - xplant[i], dy = y - yplant[i],
        marks = marks[i, ], par = effpar)
    for(i in 1:npoints(plants)) {
        phi <- outer(yrow, xcol, infi)
        if(is.infinite(asym)) part <- (phi >= denom) * (phi > 0)
        else if(asym == 1) part <- phi / (denom + (denom == 0))
        else if(asym > 0) part <- phi^asym / (denom + (denom == 0))
                # avoid division by 0
        else part <- (phi > 0) / (denom + (denom == 0)) # asym == 0
        eff <- outer(yrow, xcol, effi)
        assim <- eff * part * resmatrix
        total <- sum(assim)
        marks$aindex[i] <- total * pixarea
        if(afree) {
            marks$afree[i] <- sum(eff * (phi > 0) * resmatrix) * pixarea
        }
        if(centroid) {
            marks$cx[i] <- sum(xcol * t(assim)) / total
            marks$cy[i] <- sum(yrow * assim) / total
        }
    }
    # Add results to marks, return
    marks(plants) <- marks
    return(plants)
}
