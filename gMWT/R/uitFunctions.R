# Version: 30-11-2012, Daniel Fischer

uitPTest <- function(x,y,z,nper){
    Nx <- length(x)
    Ny <- length(y)
    Nz <- length(z)

    permTestValues <- c(rep(0,nper))

    for(i in 1:nper)
     {
	permvalues <- sample(c(x,y,z))
	permTestValues[i] <- uit.C(permvalues[1:Nx],permvalues[(Nx+1):(Nx+Ny)],permvalues[(Nx+Ny+1):(Nx+Ny+Nz)])
    }
    return(permTestValues)
}

uitATest <- function(x,y,z,obs){
  0.5 * (1-pchisq(obs,1)) + acos(getRho.C(x,y,z))/(2*pi)*(1-pchisq(obs,2))
}