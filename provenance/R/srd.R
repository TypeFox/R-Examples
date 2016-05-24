correct <- function(comp,dens,target){
    if (all(target<dens))
        stop("Target density is lower than any of the component's densities")
    len <- length(comp)
    SRD <- sum(dens*comp)/sum(comp)
    antiSRD <- sum(comp)/sum(comp/dens)
    if (target < SRD) { # if sediment SRD is higher than the target value
        return(comp*antiSRD/dens) # apply SRD correction
    } else {
        return(comp)
    }
}

restore.composition <- function(X,dens,target){
    out <- as.matrix(X,nrow=1)
    for (i in 2:30){ # run up to 30 iterations
        corrected <- correct(out[(i-1),],as.matrix(dens,nrow=1),target)
        if (sum(corrected-out[i-1,])==0){
            return(out)
        } else {
            out <- rbind(out,corrected)
        }
    }
    stop("Failed to restore the SRD composition")
}

#' Undo the effect of hydraulic sorting
#'
#' Restore the detrital composition back to a specified source rock density (SRD)
#' @param X an object of class \code{compositional}
#' @param dens a vector of rock and mineral densities
#' @param target the target density (in g/cm3)
#' @return an object of class \code{SRDcorrected}, i.e. an object of class
#' \code{compositional} which is a daughter of class \code{compositional}
#' containing the restored composition, plus one additional member called
#' \code{restoration}, containing the intermediate steps of the SRD correction
#' algorithm.
#' @examples
#' data(Namib,densities)
#' rescomp <- restore(Namib$PTHM,densities,2.71)
#' HMcomp <- c("zr","tm","rt","sph","ap","ep","gt",
#'             "st","amp","cpx","opx")
#' amcomp <- amalgamate(rescomp,Plag="P",HM=HMcomp,Opq="opaques")
#' plot(ternary(amcomp),showpath=TRUE)
#' @author Alberto Resentini and Pieter Vermeesch
#' @references Garzanti E, Ando, S and Vezzoli, G.  "Settling
#' equivalence of detrital minerals and grain-size dependence of
#' sediment composition." Earth and Planetary Science Letters 273.1
#' (2008): 138-151.
#' @seealso minsorting
#' @export
restore <- function(X,dens,target=2.71){
    if (!methods::is(X,"compositional")) stop("Input is not of class compositional")
    mydens <- get.densities(X,dens)
    out <- X
    out$restoration <- list()
    snames <- names(X)
    for (i in 1:length(snames)){
        restoration <- restore.composition(X$x[i,],mydens,target)
        out$restoration[[snames[i]]] <- restoration
        out$x[i,] <- restoration[nrow(restoration),]
    }
    class(out) <- append("SRDcorrected",class(out))
    return(out)
}
