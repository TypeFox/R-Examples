"baseCortime" <-
function (data, basev) 
{
    psisim <- data@psi.df
    tmin <- max( basev[1], 1)
    tmax <- min( basev[2], data@nt )
    lmin <- max( basev[3], 1)
    lmax <- min( basev[4], data@nl ) 
    
    for (j in tmin:tmax) {
        baseline <- sum(psisim[j, lmin:lmax])/(length(lmin:lmax))
        psisim[j, ] <- psisim[j, ] - baseline
    }
    data@psi.df <- psisim
    data
}

