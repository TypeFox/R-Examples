.packageName <- "NScluster"

SimulateIP <- function( seeds, pa, Ty=1, pmax=100, omax=3000, plot=TRUE ) {
# seeds                  #  initial seeds for a sequence of uniform random numbers
  ix <- seeds[1]
  iy <- seeds[2]
  iz <- seeds[3]
# pa                     #  the parameter values (\mu, \nu) for the random variable Poisson
  amu <- pa[1]
  anu <- pa[2]
  p <- pa[3]             #  a decay order with respect to the distance
  c <- pa[4]             #  scaling factor with respect to the distance
# Ty                     #  the variable for the standardized coordinates of points in the rectranglar region

  z <- .Call("simIP",
	as.integer(ix),
	as.integer(iy),
	as.integer(iz),
	as.double(Ty),
	as.double(amu),
	as.double(anu),
	as.double(p),
	as.double(c),
	as.integer(pmax),
	as.integer(omax))

  ier <- z[[7L]]
  if( ier == -1 ) {
    cat(sprintf(" caution : number of parents is greater than pmax %i\n", pmax))
  } else if ( ier == -2 ) {
    cat(sprintf(" caution : number of offspring is greater than omax %i\n", omax))
  } else {

  npts <- z[[1L]]
  ncl <- z[[2L]]

  parents.x <- z[[3L]][1:npts]
  parents.y <- z[[4L]][1:npts]
  parents.xy <- list(x=parents.x, y=parents.y)

  xcl <- array(z[[5L]], dim=c(pmax,omax))
  ycl <- array(z[[6L]], dim=c(pmax,omax))

  offspring.x <- NULL
  offspring.y <- NULL
  np <- 0
  for( i in 1:npts ) {
    offspring.x <- c(offspring.x, xcl[i,1:ncl[i]])
    offspring.y <- c(offspring.y, ycl[i,1:ncl[i]])
    np <- np+ncl[i]
  }
  offspring.xy <- list(x=offspring.x, y=offspring.y)

  if( plot == TRUE ) {
    old.par <- par(no.readonly=TRUE)
    par(mfrow=c(1,2), pch=19, cex=0.5, xaxs="i", yaxs="i", pty="s")

    plot(parents.x, parents.y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    title(main="Parent points",
         sub=substitute(paste("Simulated parent points  ",mu == v1), list(v1=amu)))

    plot(offspring.x, offspring.y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    title(main="Inverse-power type model", 
          sub=substitute(paste("Simulated invers-power type model  (",mu,",",nu,",p,c)" == "(",v1,",",v2,",",v3,",",v4,")"), list(v1=amu, v2=anu, v3=p, v4=c)))

    par(mfrow=old.par$mfrow, pch=old.par$pch, cex=old.par$cex, xaxs=old.par$xaxs, yaxs=old.par$yaxs, pty=old.par$pty)
  }

  SimulateIP.out <- list( n.parents=npts, parents=parents.xy, n.offspring=np, offspring=offspring.xy )
  class( SimulateIP.out ) <- "Simulate"
  return( SimulateIP.out )
  }
}


SimulateThomas <- function( seeds, pa, Ty=1, pmax=100, omax=3000, plot=TRUE )  {
# seeds                  #  initial seeds for a sequence of uniform random numbers
  ix <- seeds[1]
  iy <- seeds[2]
  iz <- seeds[3]
# pa                     #  the parameter values (\mu, \nu, \sigma) for the random variable Poisson
  amu <- pa[1]
  anu <- pa[2]
  sig <- pa[3]
# Ty                     #  the variable for the standardized coordinates of points in the rectranglar region

  z <- .Call("simThom",
	as.integer(ix),
	as.integer(iy),
	as.integer(iz),
	as.double(Ty),
	as.double(amu),
	as.double(anu),
	as.double(sig),
	as.integer(pmax),
	as.integer(omax))

  ier <- z[[7L]]
  if( ier == -1 ) {
    cat(sprintf(" caution : number of parents is greater than pmax %i\n", pmax))
  } else if ( ier == -2 ) {
    cat(sprintf(" caution : number of offspring is greater than omax %i\n", omax))
  } else {

  npts <- z[[1L]]
  ncl <- z[[2L]]
  xcl <- array(z[[5L]], dim=c(pmax,omax))
  ycl <- array(z[[6L]], dim=c(pmax,omax))

  parents.x <- z[[3L]][1:npts]
  parents.y <- z[[4L]][1:npts]
  parents.xy <- list(x=parents.x, y=parents.y)

  offspring.x <- NULL
  offspring.y <- NULL
  np <- 0
  for( i in 1:npts ) {
    offspring.x <- c(offspring.x, xcl[i,1:ncl[i]])
    offspring.y <- c(offspring.y, ycl[i,1:ncl[i]])
    np <- np+ncl[i]
  }
  offspring.xy <- list(x=offspring.x, y=offspring.y)

  if( plot == TRUE ) {
    old.par <- par(no.readonly=TRUE)
    par(mfrow=c(1,2), pch=19, cex=0.5, xaxs="i", yaxs="i", pty="s")

    plot(parents.x, parents.y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    title(main="Parent points",
          sub=substitute(paste("Simulated parent points  ",mu == v1), list(v1=amu)))

    plot(offspring.x, offspring.y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    title(main="Thomas model",
          sub=substitute(paste("Simulated Thomas model  (",mu,",",nu,",",sigma,")" == "(",v1,",",v2,",",v3,")"), list(v1=amu,v2=anu,v3=sig)))

    par(mfrow=old.par$mfrow, pch=old.par$pch, cex=old.par$cex, xaxs=old.par$xaxs, yaxs=old.par$yaxs, pty=old.par$pty)
  }

  SimulateThomas.out <- list( n.parents=npts, parents=parents.xy, n.offspring=np, offspring=offspring.xy )
  class( SimulateThomas.out ) <- "Simulate"
  return( SimulateThomas.out )
  }
}


SimulateTypeA <- function( seeds, pa, Ty=1, pmax=100, omax=3000, plot=TRUE ) {
# seeds                  #  initial seeds for a sequence of uniform random numbers
  ix <- seeds[1]
  iy <- seeds[2]
  iz <- seeds[3]
# pa                     #  the parameter values (\mu, \nu, a, \sigma1, \sigma2)
  amu <- pa[1]           #  for the random variable Poisson
  anu <- pa[2]
  a   <- pa[3]
  sig1 <- pa[4]
  sig2 <- pa[5]
# Ty                     #  the variable for the standardized coordinates of points in the rectranglar region

  z <- .Call("simA",
	as.integer(ix),
	as.integer(iy),
	as.integer(iz),
	as.double(Ty),
	as.double(amu),
	as.double(anu),
	as.double(a),
	as.double(sig1),
	as.double(sig2),
	as.integer(pmax),
	as.integer(omax))


  ier <- z[[7L]]
  if( ier == -1 ) {
    cat(sprintf(" caution : number of parents is greater than pmax %i\n", pmax))
  } else if ( ier == -2 ) {
    cat(sprintf(" caution : number of offspring is greater than omax %i\n", omax))
  } else {

    npts <- z[[1L]]
    ncl <- z[[2L]]
    xcl <- array(z[[5L]], dim=c(pmax,omax))
    ycl <- array(z[[6L]], dim=c(pmax,omax))

    parents.x <- z[[3L]][1:npts]
    parents.y <- z[[4L]][1:npts]
    parents.xy <- list(x=parents.x, y=parents.y)

    offspring.x <- NULL
    offspring.y <- NULL
    np <- 0
    for( i in 1:npts ) {
      offspring.x <- c(offspring.x, xcl[i,1:ncl[i]])
      offspring.y <- c(offspring.y, ycl[i,1:ncl[i]])
      np <- np+ncl[i]
    }
    offspring.xy <- list(x=offspring.x, y=offspring.y)

    if( plot == TRUE ) {
      old.par <- par(no.readonly=TRUE)
      par(mfrow=c(1,2), pch=19, cex=0.5, xaxs="i", yaxs="i", pty="s")

      plot(parents.x, parents.y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
      title(main="Parent points",
            sub=substitute(paste("Simulated parent points  ",mu == v1), list(v1=amu)))

      plot(offspring.x, offspring.y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
      title(main="Type A model",
            sub=substitute(paste("Simulated Type A model  (",mu,",",nu,",",a,",",sigma[1],",",sigma[2],")" == "(",v1,",",v2,",",v3,",",v4,",",v5,")"), list(v1=amu,v2=anu,v3=a,v4=sig1,v5=sig2)))

      par(mfrow=old.par$mfrow, pch=old.par$pch, cex=old.par$cex, xaxs=old.par$xaxs, yaxs=old.par$yaxs, pty=old.par$pty)
    }

    SimulateTypeA.out <- list( n.parents=npts, parents=parents.xy, n.offspring=np, offspring=offspring.xy )
    class( SimulateTypeA.out ) <- "Simulate"
    return( SimulateTypeA.out )
  }
}


SimulateTypeB <- function( seeds, pa, Ty=1, pmax=100, omax=3000, plot=TRUE ) {
# seeds                  #  initial seeds for a sequence of uniform random numbers
  ix <- seeds[1]
  iy <- seeds[2]
  iz <- seeds[3]
# pa                     #  the parameter values (\mu1, \mu2, \nu, \sigma1, \sigma2)
  amu1 <- pa[1]          #  for the random variable Poisson
  amu2 <- pa[2]
  anu  <- pa[3]
  sig1 <- pa[4]
  sig2 <- pa[5]
# Ty                     #  the variable for the standardized coordinates of points in the rectranglar region

  z <- .Call("simB",
	as.integer(ix),
	as.integer(iy),
	as.integer(iz),
	as.double(Ty),
	as.double(amu1),
	as.double(amu2),
	as.double(anu),
	as.double(sig1),
	as.double(sig2),
	as.integer(pmax),
	as.integer(omax))

  if( z[[13L]] == -1 || z[[13L]] == -2) {
    cat(sprintf(" caution : number of parents is greater than pmax %i\n", pmax))
  } else if ( z[[13L]] == -11 || z[[13L]] == -22 ) {
    cat(sprintf(" caution : number of offspring is greater than omax %i\n", omax))
  } else {

  m1 <- z[[1L]]
  m2 <- z[[7L]]
  parents1.x <- z[[3L]][1:m1]
  parents1.y <- z[[4L]][1:m1]
  parents2.x <- z[[9L]][1:m2]
  parents2.y <- z[[10L]][1:m2]
  parents.xy <- list(x=c(parents1.x, parents2.x), y=c(parents1.y, parents2.y))

  np1 <- sum(z[[2L]][1:m1])
  offspring1.x <- z[[5L]][1:np1]
  offspring1.y <- z[[6L]][1:np1]
  np2 <- sum(z[[8L]][1:m2])
  offspring2.x <- z[[11L]][1:np2]
  offspring2.y <- z[[12L]][1:np2]
  offspring.xy <- list(x=c(offspring1.x, offspring2.x), y=c(offspring1.y, offspring2.y))

  if( plot == TRUE ) {
    old.par <- par(no.readonly=TRUE)
    par(mfrow=c(1,2), pch=19, cex=0.5, xaxs="i", yaxs="i", pty="s")

    x <- parents.xy$x
    y <- parents.xy$y
    plot(x, y, type="p", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    title(main="Parent points",
          sub=substitute(paste("Simulated parent points  ",mu[1] == v1, " and ", mu[2] == v2), list(v1=amu1,v2=amu2)))

    x <- offspring.xy$x
    y <- offspring.xy$y
    plot(x, y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    title(main="Type B model",
          sub=substitute(paste("Simulsted TypeB model  (",mu[1],",",mu[2],",",nu,",",sigma[1],",",sigma[2],")" == "(",v1,",",v2,",",v3,",",v4,",",v5,")"), list(v1=amu1,v2=amu2,v3=anu,v4=sig1,v5=sig2)))

    par(mfrow=old.par$mfrow, pch=old.par$pch, cex=old.par$cex, xaxs=old.par$xaxs, yaxs=old.par$yaxs, pty=old.par$pty)
  }

  SimulateTypeB.out <- list( n.parents1=m1, n.offspring1=np1, n.parents2=m2, n.offspring2=np2, parents=parents.xy, offspring=offspring.xy )
  class( SimulateTypeB.out ) <- "Simulate2"
  return( SimulateTypeB.out )
  }
}


SimulateTypeC <- function( seeds, pa1, pa2, Ty=1, pmax=100, omax=3000, plot=TRUE ) {
# seeds                  #  initial seeds for a sequence of uniform random numbers
  ix <- seeds[1]
  iy <- seeds[2]
  iz <- seeds[3]
# pa1                    #  the parameter values (\mu_1, \nu_1, \sigma_1)
  amu1 <- pa1[1]         #  for the random variable Poisson
  anu1 <- pa1[2]
  sig1 <- pa1[3]
# pa2                    #  the parameter values (\mu_2, \nu_2, \sigma_2)
  amu2 <- pa2[1]         #  for the random variable Poisson
  anu2 <- pa2[2]
  sig2 <- pa2[3]
# Ty                     #  the variable for the standardized coordinates of points in the rectranglar region

  z <- .Call("simC",
	as.integer(ix),
	as.integer(iy),
	as.integer(iz),
	as.double(Ty),
	as.double(amu1),
	as.double(amu2),
	as.double(anu1),
	as.double(anu2),
	as.double(sig1),
	as.double(sig2),
	as.integer(pmax),
	as.integer(omax))

  if( z[[13L]] == -1 || z[[13L]] == -2) {
    cat(sprintf(" caution : number of parents is greater than pmax %i\n", pmax))
  } else if ( z[[13L]] == -11 || z[[13L]] == -22 ) {
    cat(sprintf(" caution : number of offspring is greater than omax %i\n", omax))
  } else {

  m1 <- z[[1L]]
  m2 <- z[[7L]]
  parents1.x <- z[[3L]][1:m1]
  parents1.y <- z[[4L]][1:m1]
  parents2.x <- z[[9L]][1:m2]
  parents2.y <- z[[10L]][1:m2]
  parents.xy <- list(x=c(parents1.x, parents2.x), y=c(parents1.y, parents2.y))

  np1 <- sum(z[[2L]][1:m1])
  offspring1.x <- z[[5L]][1:np1]
  offspring1.y <- z[[6L]][1:np1]
  np2 <- sum(z[[8L]][1:m2])
  offspring2.x <- z[[11L]][1:np2]
  offspring2.y <- z[[12L]][1:np2]
  offspring.xy <- list(x=c(offspring1.x, offspring2.x), y=c(offspring1.y, offspring2.y))

  if( plot == TRUE ) {
    old.par <- par(no.readonly=TRUE)
    par(mfrow=c(1,2), pch=19, cex=0.5, xaxs="i", yaxs="i", pty="s")

    x <- parents.xy$x
    y <- parents.xy$y
    plot(x, y, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    title(main="Parent points",
          sub=substitute(paste("Simulated parent points  ",mu[1] == v1, " and ", mu[2] == v2), list(v1=amu1,v2=amu2)))

    x <- offspring.xy$x
    y <- offspring.xy$y
    plot(x, y, xlim=c(0,1), ylim=c(0,1),xlab="", ylab="")
    title(main="Type C model",
          sub=substitute(paste("Simulated Type C model  (",mu[1],",",mu[2],",",nu[1],",",nu[2],",",sigma[1],",",sigma[2],")" == "(",v1,",",v2,",",v3,",",v4,",",v5,",",v6,")"), list(v1=amu1,v2=amu2,v3=anu1,v4=anu2,v5=sig1,v6=sig2)))

    par(mfrow=old.par$mfrow, pch=old.par$pch, cex=old.par$cex, xaxs=old.par$xaxs, yaxs=old.par$yaxs, pty=old.par$pty)
  }

  SimulateTypeC.out <- list( n.parents1=m1, n.offspring1=np1, n.parents2=m2, n.offspring2=np2, parents=parents.xy, offspring=offspring.xy )
  class( SimulateTypeC.out ) <- "Simulate2"
  return( SimulateTypeC.out )
  }
}


SimplexThomas <- function( offspring, pa, Ty=1, eps=0.1e-2, process=0, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu, \nu, \sigma) for the random variable Poisson
  sclmu <- pa[1]
  sclnu <- pa[2]
  scls <- pa[3]
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

# parameter
  ipflg <- process
  if( process==0 && plot==TRUE ) ipflg <- 2
  if( process==1 && plot==TRUE ) ipflg <- 3

  n <- 3
  itmax <- 1000
  itmax1 <- 1
  if ( ipflg > 1 ) itmax1 <- itmax+1
  ipmax <- itmax*2
  if( (ipflg==0) || (ipflg==2) ) ipmax <- 1

  z <- .Call("smplxThom",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(ty),
	as.double(sclmu),
	as.double(sclnu),
	as.double(scls),
	as.double(eps),
	as.integer(itmax),
	as.integer(itmax1),
	as.integer(ipmax),
	as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array( z[[2L]], dim=c(ipmax,n))
  if( nip == 1 ) mples <- data.frame(mples) 
  if( nip != 1 ) mples <- data.frame(mples[1:nip,1:n])
  names(mples) <- c("mu", "nu", "sigma") 

  it1 <- z[[6L]] + 1
  it2 <- 1
  if ( (ipflg==2) || (ipflg==3) ) it2 <- it1
  xx <- array(z[[3L]], dim=c(n,itmax1))

  if ( process < 2 ) {
    f <- z[[5L]][it2]
    para <- xx[1:n,it2]
    std <- z[[4L]][it2]
  } else {
    f <- z[[5L]][1:it2]
    para <- xx[1:n,1:it2]
    std <- z[[4L]][1:it2]
  }

  if( plot == TRUE ) {
    param <- xx[1:n,1:it2]
    stddev <- z[[4L]][1:it2]
    plot.Simplex("T", stddev, param)
  }

  if( (process==0 ) || (process==2) )  SimplexThomas.out <- list( pa.init=pa, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para )

  if( (process==1) || (process==3) )  SimplexThomas.out <- list( pa.init=pa, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para, ipri=ipri )

  class( SimplexThomas.out ) <- "Simplex"
  return( SimplexThomas.out )
}


SimplexIP <- function( offspring, pa, Ty=1, x2, skip=1, eps=0.1e-2, process=0, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu, \nu, \sigma) for the random variable Poisson
  sclmu <- pa[1]
  sclnu <- pa[2]
  sclp <- pa[3]
  sclc <- pa[4]
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region
# x2
# skip                   

# parameter
  ipflg <- process
  if( process==0 && plot==TRUE ) ipflg <- 2
  if( process==1 && plot==TRUE ) ipflg <- 3

  n <- 4
  itmax <- 1000
  itmax1 <- 1
  if ( ipflg > 1 ) itmax1 <- itmax+1
  ipmax <- itmax*2
  if( (ipflg==0) || (ipflg==2) ) ipmax <- 1

  z <- .Call("smplxIP",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.integer(skip),
	as.double(ty),
	as.double(sclmu),
	as.double(sclnu),
	as.double(sclp),
	as.double(sclc),
	as.double(x2),
	as.double(eps),
	as.integer(itmax),
	as.integer(itmax1),
	as.integer(ipmax),
	as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array( z[[2L]], dim=c(ipmax,n))
  if( nip == 1 ) mples <- data.frame(mples) 
  if( nip != 1 ) mples <- data.frame(mples[1:nip,1:n])
  names(mples) <- c("mu", "nu", "a", "p") 

  it1 <- z[[6L]] + 1
  it2 <- 1
  if ( (ipflg==2) || (ipflg==3) ) it2 <- it1
  xx <- array(z[[3L]], dim=c(n,itmax1))

  if ( process < 2 ) {
    f <- z[[5L]][it2]
    para <- xx[1:n,it2]
    std <- z[[4L]][it2]
  } else {
    f <- z[[5L]][1:it2]
    para <- xx[1:n,1:it2]
    std <- z[[4L]][1:it2]
  }

  if( plot == TRUE ) {
    param <- xx[1:n,1:it2]
    stddev <- z[[4L]][1:it2]
    plot.Simplex("IP", stddev, param)
  }

  if( (process==0 ) || (process==2) )  SimplexIP.out <- list( pa.init=pa, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para )

  if( (process==1) || (process==3) )  SimplexIP.out <- list( pa.init=pa, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para, ipri=ipri )

  class( SimplexIP.out ) <- "Simplex"
  return( SimplexIP.out )
}


SimplexTypeA <- function( offspring, pa, Ty=1, x2, skip=1, eps=0.1e-2, process=0, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu, \nu, \sigma) for the random variable Poisson
  sclmu <- pa[1]
  sclnu <- pa[2]
  scla <- pa[3]
  scls1 <- pa[4]
  scls2 <- pa[5]
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region
# x2
# skip                   

# parameter
  ipflg <- process
  if( process==0 && plot==TRUE ) ipflg <- 2
  if( process==1 && plot==TRUE ) ipflg <- 3

  n <- 5
  itmax <- 1000
  itmax1 <- 1
  if ( ipflg > 1 ) itmax1 <- itmax+1
  ipmax <- itmax*2
  if( (ipflg==0) || (ipflg==2) ) ipmax <- 1

  z <- .Call("smplxA", 
	as.double(x),
	as.double(y),
	as.integer(np),
	as.integer(skip),
	as.double(ty),
	as.double(sclmu),
	as.double(sclnu),
	as.double(scla),
	as.double(scls1),
	as.double(scls2),
	as.double(x2),
	as.double(eps),
	as.integer(itmax),
	as.integer(itmax1),
	as.integer(ipmax),
	as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array( z[[2L]], dim=c(ipmax,n))
  if( nip == 1 ) mples <- data.frame(mples) 
  if( nip != 1 ) mples <- data.frame(mples[1:nip,1:n])
  names(mples) <- c("mu", "nu", "a", "sigma1", "sigma2") 

  it1 <- z[[6L]] + 1
  it2 <- 1
  if ( (ipflg==2) || (ipflg==3) ) it2 <- it1
  xx <- array(z[[3L]], dim=c(n,itmax1))

  if ( process < 2 ) {
    f <- z[[5L]][it2]
    para <- xx[1:n,it2]
    std <- z[[4L]][it2]
  } else {
    f <- z[[5L]][1:it2]
   para <- xx[1:n,1:it2]
    std <- z[[4L]][1:it2]
  }

  if( plot == TRUE ) {
    param <- xx[1:n,1:it2]
    stddev <- z[[4L]][1:it2]
    plot.Simplex("A", stddev, param)
  }

  if( (process==0 ) || (process==2) )  SimplexTypeA.out <- list( pa.init=pa, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para )

  if( (process==1) || (process==3) )  SimplexTypeA.out <- list( pa.init=pa, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para, ipri=ipri )

  class( SimplexTypeA.out ) <- "Simplex"
  return( SimplexTypeA.out )
}


SimplexTypeB <- function( offspring, pa, Ty=1, eps=0.1e-2, process=0, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu, \nu, \sigma) for the random variable Poisson
  sclamu1 <- pa[1]
  sclamu2 <- pa[2]
  sclanu <- pa[3]
  scls1 <- pa[4]
  scls2 <- pa[5]
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

# parameter
  ipflg <- process
  if( process==0 && plot==TRUE ) ipflg <- 2
  if( process==1 && plot==TRUE ) ipflg <- 3

  n <- 5
  itmax <- 1000
  itmax1 <- 1
  if ( ipflg > 1 ) itmax1 <- itmax+1
  ipmax <- itmax*2
  if( (ipflg==0) || (ipflg==2) ) ipmax <- 1

  z <- .Call("smplxB",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(ty),
	as.double(sclamu1),
	as.double(sclamu2),
	as.double(sclanu),
	as.double(scls1),
	as.double(scls2),
	as.double(eps),
	as.integer(itmax),
	as.integer(itmax1),
	as.integer(ipmax),
	as.integer(ipflg))

  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array( z[[2L]], dim=c(ipmax,n))
  if( nip == 1 ) mples <- data.frame(mples) 
  if( nip != 1 ) mples <- data.frame(mples[1:nip,1:n])
  names(mples) <- c("mu", "nu", "a", "sigma1", "sigma2") 

  it1 <- z[[6L]] + 1
  it2 <- 1
  if ( (ipflg==2) || (ipflg==3) ) it2 <- it1
  xx <- array(z[[3L]], dim=c(n,itmax1))

  if ( process < 2 ) {
    f <- z[[5L]][it2]
    para <- xx[1:n,it2]
    std <- z[[4L]][it2]
  } else {
    f <- z[[5L]][1:it2]
    para <- xx[1:n,1:it2]
    std <- z[[4L]][1:it2]
  }

  if( plot == TRUE ) {
    param <- xx[1:n,1:it2]
    stddev <- z[[4L]][1:it2]
    plot.Simplex("B", stddev, param)
  }

  pini <- c(pa[1]+pa[2], pa[3], pa[1]/(pa[1]+pa[2]), pa[4], pa[5])

  if( (process==0 ) || (process==2) )  SimplexTypeB.out <- list( pa.init=pini, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para )

  if( (process==1) || (process==3) )  SimplexTypeB.out <- list( pa.init=pini, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para, ipri=ipri )

  class( SimplexTypeB.out ) <- "Simplex"
  return( SimplexTypeB.out )
}


SimplexTypeC <- function( offspring, pa1, pa2, Ty=1, eps=0.1e-2, process=0, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu_1,\mu_2,\nu_1,\nu_2,\sigma_1,sigma_2) for the random variable Poisson
  amu1 <- pa1[1]
  amu2 <- pa2[1]
  anu1 <- pa1[2]
  anu2 <- pa2[2]
  scls1 <- pa1[3]
  scls2 <- pa2[3]
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

# parameter
  ipflg <- process
  if( process==0 && plot==TRUE ) ipflg <- 2
  if( process==1 && plot==TRUE ) ipflg <- 3

  n <- 5
  itmax <- 1000
  itmax1 <- 1
  if ( ipflg > 1 ) itmax1 <- itmax+1
  ipmax <- itmax*2
  if( (ipflg==0) || (ipflg==2) ) ipmax <- 1

 z <- .Call("smplxC",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(ty),
	as.double(amu1),
	as.double(amu2),
	as.double(anu1),
	as.double(anu2),
	as.double(scls1),
	as.double(scls2),
	as.double(eps),
	as.integer(itmax),
	as.integer(itmax1),
	as.integer(ipmax),
	as.integer(ipflg))

#  mple <- alam, anu, a, sigma1, sigma2
  nip <- z[[7L]]
  ipri <- z[[8L]][1:nip]
  mples <- array( z[[2L]], dim=c(ipmax,n))
  if( nip == 1 ) mples <- data.frame(mples) 
  if( nip != 1 ) mples <- data.frame(mples[1:nip,1:n])
  names(mples) <- c("lambda", "nu", "a", "sigma1", "sigma2") 

  it1 <- z[[6L]] + 1
  it2 <- 1
  if ( (ipflg==2) || (ipflg==3) ) it2 <- it1
  xx <- array(z[[3L]], dim=c(n,itmax1))

  if ( process < 2 ) {
    f <- z[[5L]][it2]
    para <- xx[1:n,it2]
    std <- z[[4L]][it2]
  } else {
    f <- z[[5L]][1:it2]
    para <- xx[1:n,1:it2]
    std <- z[[4L]][1:it2]
  }

  if( plot == TRUE ) {
    param <- xx[1:n,1:it2]
    stddev <- z[[4L]][1:it2]
    plot.Simplex("C", stddev, param)
  }

  scllam <- amu1*anu1+amu2*anu2
  scla <- amu1*anu1/scllam
  pini <- c(scllam, anu1, scla, pa1[3], pa2[3])

  if( (process==0 ) || (process==2) )  SimplexTypeC.out <- list( pa.init=pini, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para )

  if( (process==1) || (process==3) )  SimplexTypeC.out <- list( pa.init=pini, logL.p=z[[1L]][1:nip], mple=mples, logL.s=f, stderr=std, pa.normal=para, ipri=ipri )

  class( SimplexTypeC.out ) <- "Simplex"
  return( SimplexTypeC.out )
}


PalmThomas <- function( offspring, pa, delta, Ty, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu_i,\nu_i,\sigma_i)
  m <- dim(pa)[1]
  amu <- rep(0,m); anu <- rep(0,m); v <- rep(0,m)
  for( i in 1:m ) {
    amu[i] <- pa[i,1]
    anu[i] <- pa[i,2]
    v[i] <- pa[i,3]
  }
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

  R=0.5e0
  rmax=R/delta
  jmax <- as.integer(rmax)

  z <- .Call("palmT",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(delta),
	as.double(ty),
	as.double(amu),
	as.double(anu),
	as.double(v),
	as.integer(m),
	as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax)*delta
  palm1 <- array(z[[2L]], dim=c(jmax,m))

  PalmThomas.out <- list( r=r, np.palm=palm, palm.normal=palm1 )

  if( plot == TRUE )  {
    plot.Palm("T", r, palm, palm1)
    return( invisible(PalmThomas.out) )
  } else {
    class( PalmThomas.out ) <- "Palm"
    return( PalmThomas.out )
  }
}


PalmIP <- function( offspring, pa, delta, Ty, x2, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu_i,\nu_i,p_i,c_i)
  m <- dim(pa)[1]
  amu <- rep(0,m); anu <- rep(0,m); p <- rep(0,m); c <- rep(0,m)
  for( i in 1:m ) {
    amu[i] <- pa[i,1]
    anu[i] <- pa[i,2]
    p[i] <- pa[i,3]
    c[i] <- pa[i,4]
  }
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

  R=0.5e0
  rmax=R/delta
  jmax <- as.integer(rmax)

  z <- .Call("palmIP",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(delta),
	as.double(ty),
	as.double(x2),
	as.double(amu),
	as.double(anu),
	as.double(p),
	as.double(c),
	as.integer(m),
	as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax)*delta
  palm1 <- array(z[[2L]], dim=c(jmax,m))

  PalmIP.out <- list( r=r, np.palm=palm, palm.normal=palm1 )

  if( plot == TRUE )  {
    plot.Palm("IP", r, palm, palm1)
    return( invisible(PalmIP.out) )
  } else {
    class( PalmIP.out ) <- "Palm"
    return( PalmIP.out )
  }
}


PalmTypeA <- function( offspring, pa, delta, Ty, x2, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu_i,\nu_i,a,\sigma_{1,i},\sigma_{2,i})
  m <- dim(pa)[1]
  amu <- rep(0,m); anu <- rep(0,m); a <- rep(0,m); s1 <- rep(0,m); s2 <- rep(0,m)
  for( i in 1:m ) {
    amu[i] <- pa[i,1]
    anu[i] <- pa[i,2]
    a[i] <- pa[i,3]
    s1[i] <- pa[i,4]
    s2[i] <- pa[i,5]
  }
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

  R=0.5e0
  rmax=R/delta
  jmax <- as.integer(rmax)

  z <- .Call("palmA",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(delta),
	as.double(ty),
	as.double(x2),
	as.double(amu),
	as.double(anu),
	as.double(a),
	as.double(s1),
	as.double(s2),
	as.integer(m),
	as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax)*delta
  palm1 <- array(z[[2L]], dim=c(jmax,m))

  PalmTypeA.out <- list( r=r, np.palm=palm, palm.normal=palm1 )

  if( plot == TRUE )  {
    plot.Palm("A", r, palm, palm1)
    return( invisible(PalmTypeA.out) )
  } else {
    class( PalmTypeA.out ) <- "Palm"
    return( PalmTypeA.out )
  }
}


PalmTypeB <- function( offspring, pa, delta, Ty, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\mu_i,\nu_i,a, \sigma_{1,i}, \sigma_{2,i})
  m <- dim(pa)[1]
  amu <- rep(0,m); anu <- rep(0,m); a <- rep(0,m); s1 <- rep(0,m); s2 <- rep(0,m)
  for( i in 1:m ) {
    amu[i] <- pa[i,1]
    anu[i] <- pa[i,2]
    a[i] <- pa[i,3]
    s1[i] <- pa[i,4]
    s2[i] <- pa[i,5]
  }
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

  R=0.5e0
  rmax=R/delta
  jmax <- as.integer(rmax)

  z <- .Call("palmB",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(delta),
	as.double(ty),
	as.double(amu),
	as.double(anu),
	as.double(a),
	as.double(s1),
	as.double(s2),
	as.integer(m),
	as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax)*delta
  palm1 <- array(z[[2L]], dim=c(jmax,m))

  PalmTypeB.out <- list( r=r, np.palm=palm, palm.normal=palm1 )

  if( plot == TRUE )  {
    plot.Palm("B", r, palm, palm1)
    return( invisible(PalmTypeB.out) )
  } else {
    class( PalmTypeB.out ) <- "Palm"
    return( PalmTypeB.out )
  }
}


PalmTypeC <- function( offspring, pa, delta, Ty, plot=TRUE ) {

  x <- offspring$x
  y <- offspring$y
  np <- length(x)
# pa                     #  the parameter values (\lambda_i,\nu_{1,i},a, \sigma_{1,i}, \sigma_{2,i})
  m <- dim(pa)[1]
  alam <- rep(0,m); anu1 <- rep(0,m); a <- rep(0,m); s1 <- rep(0,m); s2 <- rep(0,m)
  for( i in 1:m ) {
    alam[i] <- pa[i,1]
    anu1[i] <- pa[i,2]
    a[i] <- pa[i,3]
    s1[i] <- pa[i,4]
    s2[i] <- pa[i,5]
  }
  ty <- Ty               #  the variable for the standardized coordinates of points in the rectranglar region

  R=0.5e0
  rmax=R/delta
  jmax <- as.integer(rmax)

  z <- .Call("palmC",
	as.double(x),
	as.double(y),
	as.integer(np),
	as.double(delta),
	as.double(ty),
	as.double(alam),
	as.double(anu1),
	as.double(a),
	as.double(s1),
	as.double(s2),
	as.integer(m),
	as.integer(jmax))

  palm <- z[[1L]]
  r <- rep(1:jmax)*delta
  palm1 <- array(z[[2L]], dim=c(jmax,m))

  PalmTypeC.out <- list( r=r, np.palm=palm, palm.normal=palm1 )

  if( plot == TRUE )  {
    plot.Palm("C", r, palm, palm1)
    return( invisible(PalmTypeC.out) )
  } else {
    class( PalmTypeC.out ) <- "Palm"
    return( PalmTypeC.out )
  }
}


print.Simulate <- function(x,...)
{
  cat(sprintf("\nNumber of parents =\t %i\n", x$n.parents))
  cat(sprintf("Total number of offspring =\t %i\n", x$n.offspring))
}

print.Simulate2 <- function(x,...)
{
  total.p <- x$n.parents1 + x$n.parents2
  total.o <- x$n.offspring1 + x$n.offspring2
  cat(sprintf("\nNumber of parents =\t %i\t%i\tTotal\t%i\n", x$n.parents1, x$n.parents2, total.p))
  cat(sprintf("Number of offspring =\t %i\t%i\tTotal\t%i\n", x$n.offspring1, x$n.offspring2, total.o))
}


plot.Palm <- function(type,r,palm,palm1)
{

    old.par <- par(no.readonly=TRUE)
    ymax <- max(palm, palm1)*2
    plot(x=r, y=palm, pch=20, ylim=c(0.5,ymax), log="xy", xlab="r", ylab="Palm intensity")
    if( type == "T" ) title(main = "Thomas model", sub = " dots : non-parametric / red : true / green : MPLE")
    if( type == "IP" ) title(main = "Inverse-power type model", sub = " dots : non-parametric / red : true / green : MPLE")
    if( type == "A" ) title(main = "Type A model", sub = " dots : non-parametric / red : true / green : MPLE")
    if( type == "B" ) title(main = "Type B model", sub = " dots : non-parametric / red : true / green : MPLE")
    if( type == "C" ) title(main = "Type C model", sub = " dots : non-parametric / red : true / green : MPLE")

    m = dim(palm1)[2]
    for( j in 1:m ) {
      par(new=TRUE)
      plot(x=r, y=palm1[,j], type="l", ylim=c(0.5,ymax), col=1+j, log="xy", xlab="", ylab="")
    }
    abline(h=1)

    par <- old.par
}


print.Palm <- function(x,...)
{
  cat("\n\t\t Palm intensity function\n")
  cat("\tr\t non-parametric\t     normalized Palm\n")
  n <- length(x$r)
  m <- dim(x$palm.normal)[2]
  for( i in 1:n ) {
    cat(sprintf("\n %10.3f\t%12.5f\t",x$r[i],x$np.palm[i]))
    for( j in 1:m ) cat(sprintf("\t%12.5f", x$palm.normal[i,j]))
  }
  cat("\n\n")

}

print.Simplex <- function(x,...)
{
  it2 <- 1
  if( is.matrix(x$pa.normal) == TRUE ) it2 <- dim(x$pa.normal)[2]
  n <- length(x$pa.init)
#  ipri2 <- 1
#  if( is.matrix(x$mple) == TRUE ) ipri2 <- length(x$ipri)
  ipri2 <- max(length(x$ipri), 1)

  pname <- names(x$mple)
  if( ipri2 > 1 ) {
    cat("\n\n Minimizing the negative Palm log-likelihood function\n")
    cat("\t\t-log L\t")
    for( i in 1:n ) cat(sprintf("\t%12s", pname[i]))
    cat("\n") 
    for( j in 1:ipri2 ) {
      if( x$ipri[j] == -1 ) cat(sprintf(" testfn = %18.10e", x$logL.p[j]))
      if( x$ipri[j] == 1 )  cat(sprintf(" update = %18.10e", x$logL.p[j]))
      for( i in 1:n )  cat(sprintf("\t%12.4e", x$mple[j,i]))
      cat("\n")
    }
  }

  if( it2 > 1 ) {
    cat("\n\n Optimization procedure by the simplex with the normalized parameters\n")
    cat(" stage\t\t -logL\t\t standard error")
    for( i in 1:n ) cat(sprintf("\t\t x(%i)", i))
    cat("\n")
    for( j in 1:it2 ) {
      cat(sprintf(" %i\t %15.3f\t %12.4e", j-1, x$logL.s[j], x$stderr[j]))
      for( i in 1:n ) cat(sprintf("\t %12.4e", x$pa.normal[i,j]))
      cat("\n")
    }
  }

  if( ipri2 == 1 ) mples <- x$mple
  if( ipri2 != 1 ) mples <- x$mple[ipri2,]
  if( it2 == 1 ) pa.norm <- x$pa.normal
  if( it2 != 1 ) pa.norm <- x$pa.normal[,it2]
  param <- rep(0,n)
  for( i in 1:n ) param[i] <- x$pa.init[i] * pa.norm[i]

  cat("\n\n Parameter")
  for( i in 1:n ) cat(sprintf("\t%12s", pname[i]))
  cat("\n\n") 
  cat(" True values")
  for( i in 1:n ) cat(sprintf("\t%12.4f", x$pa.init[i]))
  cat("\n\n")
  cat(" MPLE : Palm")
  for( i in 1:n )   cat(sprintf("\t%12.4f", mples[i]))
  cat("\n\n")
  cat("      : simplex")
  for( i in 1:n ) cat(sprintf("\t%12.4f", param[i]))
  cat("\n")
  cat("     normalized")
  for( i in 1:n ) cat(sprintf("\t (%10.4f)", pa.norm[i]))
  cat("\n")

}

plot.Simplex <- function(type, std, para)
{
    old.par <- par(no.readonly=TRUE)

    n <- dim(para)[1]
    if( type != "T" ) plot(std, type="l",ylim=c(0,2), sub=substitute(paste("standard deviations : black\nnormalized parameters x(1),...,x(",v1,") : color"), list(v1=n)), xlab="", ylab="")
    if( type == "T" ) plot(std, type="l",ylim=c(0,2), sub="standard deviations : black\nnormalized parameters x(1),x(2),x(3) : color", xlab="", ylab="")

    if( type == "T" ) title(main = "Thomas model")
    if( type == "IP" ) title(main = "Inverse-power type model")
    if( type == "A" ) title(main = "Type A model")
    if( type == "B" ) title(main = "Type B model")
    if( type == "C" ) title(main = "Type C model")

    for(i in 1:n) {
      par(new=TRUE)
      plot(para[i,], type="l",ylim=c(0,2), ylab="", xlab="", col=i+1)
    }

    par <- old.par
}

.noGenerics <- TRUE

options(warn.FPU=FALSE)

