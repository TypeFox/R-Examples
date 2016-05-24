# get size shift in phi units
getSS <- function(dens,medium,SRD,phi){
    ncat <- length(dens)
    SS <- matrix(,nrow=1, ncol=ncat) # vector to store size-shifts
    microns <- 2^(-phi)*1000	# transform grainsize from phi units to microns
    if (medium == "air") {
        v <- (2/3*981*(SRD-0.0012)*(microns/10000)/0.0012)^(1/2) # settling velocity
        SS <- log2((dens-0.0012)/(SRD-0.0012)) # size-shifts
    } else { # if sediment was deposited in water
        if (medium == "seawater") {
            fd<-1.025  # fluid density
            fv<-0.0105 # viscosity
        }
        if (medium == "freshwater"){
            fd<-1
            fv<-0.01
        }
        if (phi>3.5) { # if sediment is very fine --> use Stokes
            v<-981*(SRD-fd)*microns^2/18/fv
            SS<-log2((dens-fd)/(SRD-fd))/2
        } else { # if sediment is sand-sized --> use Cheng
            v<-((25+1.2*((981*(SRD-fd)*(microns/10000)^3/fv^2)^(2/3)))^(1/2)-5)^(3/2)*
                fv/(microns/10000)
            SS<-log2((dens-fd)/(SRD-fd))-3/2*
                log2((v/fv+((v/fv)^2+48*(981*(dens-fd)/fv^2)^(2/3))^(1/2))/
                     (v/fv+((v/fv)^2+48*(981*(SRD-fd)/fv^2)^(2/3))^(1/2)))	
        }
    }
    return(SS)
}

#' Assess settling equivalence of detrital components
#'
#' Models grain size distribution of minerals and rock fragments of different densities
#' @param X an object of class \code{compositional}
#' @param dens a vector of mineral and rock densities
#' @param sname sample name if unspecified, the first sample of the dataset will be used
#' @param phi the mean grain size of the sample in Krumbein's phi units
#' @param sigmaphi the standard deviation of the grain size distirbution, in phi units
#' @param medium the transport medium, one of either "air", "freshwater" or "seawater"
#' @param from the minimum grain size to be evaluated, in phi units
#' @param to the maximum grain size to be evaluated, in phi units
#' @param by the grain size interval of the output table, in phi units
#' @return an object of class \code{minsorting}, i.e. a list with two tables:
#'
#' mfract: the grain size distribution of each mineral (sum of the columns = 1)
#' 
#' mcomp: the composition of each mineral (sum of the rows = 1)
#' @author Alberto Resentini and Pieter Vermeesch
#' @references Resentini, A, Malusa, M G and Garzanti, E. "MinSORTING:
#' An Excel worksheet for modelling mineral grain-size distribution in
#' sediments, with application to detrital geochronology and
#' provenance studies." Computers & Geosciences 59 (2013): 90-97.
#'
#' Garzanti, E, Ando, S and Vezzoli, G. "Settling equivalence of
#' detrital minerals and grain-size dependence of sediment
#' composition." Earth and Planetary Science Letters 273.1 (2008):
#' 138-151.
#' @seealso restore
#' @examples
#' data(endmembers,densities)
#' distribution <- minsorting(endmembers,densities,sname='ophiolite',phi=2,
#'                            sigmaphi=1,medium="seawater",by=0.05)
#' plot(distribution,cumulative=FALSE)
#' @export
minsorting <- function(X,dens,sname=NULL,phi=2,sigmaphi=1,medium="freshwater",from=-2.25,to=5.5,by=0.25){
    if (!methods::is(X,"compositional")) stop("Input does not have class compositional")
    if (is.null(sname)) sname <- names(X)[1]
    Y <- subset(X,select=sname)
    mydens <- get.densities(Y,dens)
    ncat <- length(Y$x)
    SRD <- sum(mydens*Y$x)/sum(Y$x) # calculate sediment SRD

    SS <- getSS(mydens,medium,SRD,phi) # Size Shift
    
    M <- phi+SS # calculate mean size for each component
    w <- Y$x*mydens/SRD # calculate weight distribution
    g <- t(seq(from=from, to=to, by=by)) # grain sizes to be considered
    n <- matrix(,nrow=length(g), ncol=ncat) # matrix to store distribution
    colnames(n) <- colnames(Y$x)
    rownames(n) <- g
    
    for (i in 1:ncat) {
        n[,i] <- w[,i]*stats::pnorm(g, M[,i], sigmaphi) # pnorm=cumulative; #dnorm=gaussian
    }
    N <- diff(n) # weight amount in each class
    mfract <- t(t(N)/colSums(N))*100 # store % of each component in each size class

    r <- matrix(,nrow=length(g), ncol=ncat) # matrix to store the normal distribution
    for (i in 1:ncat) r[,i] <- w[,i]*stats::dnorm(g, M[,i], sigmaphi) # normal distribution
    y <- rowSums(r) # sum of different composition
    W <- r/y*100  # weight composition for the different size-classes
    srd <- rowSums(W %*% diag(mydens))/100 # from weight to volume compositions
    wt2vol1 <- W * srd # first step to transform from weight to volume
    wt2vol2 <- wt2vol1 %*% diag(1/mydens) # second step
    volsum <- rowSums(wt2vol2) # closure constant for each row
    mcomp <- wt2vol2/volsum*100 # model composition of each fraction
    dimnames(mcomp) <- list(g,names(Y$x))

    out <- list(mfract=mfract,mcomp=mcomp)
    class(out) <- "minsorting"
    return(out)
}
