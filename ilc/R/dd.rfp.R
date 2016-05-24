dd.rfp <-
function(ddata, rfp){
    M <- insp.dd(ddata)
    E <- insp.dd(ddata, 'pop')
    D <- M*E
    k <- dim(M)[1]; n <- dim(M)[2]
    
    # age-specific standard deviations to be applied to the exposure matrix
    # arbitrarily chosen as sqrt(SE)
    asd <- apply(E, 1, function(x) sqrt(sd(x,na.rm=T)/sqrt(length(x))))
    
    if (is.null(names(rfp))) names(rfp) <- letters[as.seq(rfp)]
    
    rM <- M
    rE <- E
    for(i in as.seq(rfp)){
        # apply a poisson error to the mortality rate matrix:
        iM <- M * matrix(rpois(k*n, exp(rfp[i])), k, n)
        # apply an age-specific normal error to the exposure matrix:
        iE <- E + matrix(rnorm(k*n, sd=asd), k, n, byrow=F)
        # reset negative and NA exposures:
        iE[bool(iE < 0, na=T)] <- 0
        rM <- cbind(rM, iM)
        rE <- cbind(rE, iE)
        D <- cbind(D, iM*iE)
    }
    M <- array(rM, dim=c(k,n,length(rfp)+1))
    E <- array(rE, dim=c(k,n,length(rfp)+1))
    D <- array(D, dim=c(k,n,length(rfp)+1))
    ret <- list(age=ddata$age, year=ddata$year, covariates = list(c('base', names(rfp))),
                deaths=D, pop=E, mu=M, type = ddata$type, label = ddata$label, name=names(ddata$rate))
    names(ret$covariates) <- 'X'
    return(structure(ret, class='rhdata'))
}
