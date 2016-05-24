
##==============================================================================
## Diffusive transport in polar coordinates (r, theta)
##==============================================================================

tran.polar <- function(C, C.r.up=NULL,    C.r.down=NULL, 
                        C.theta.up=NULL,  C.theta.down=NULL, 
                        
                        flux.r.up=NULL,   flux.r.down=NULL, 
                        flux.theta.up=NULL, flux.theta.down=NULL, 
                        
                        cyclicBnd = NULL,
                                                
                        D.r=1, D.theta=D.r,

                        r=NULL, theta=NULL,
                        
                        full.output = FALSE )
											
{
# ------------------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------------------

# a function to check the dimensionality of the system
  checkDim <- function (dr, dtheta, nx="dr", ntheta="dtheta") {
    if (length(dr) != 1 )
      stop (nx,", should have length 1")
    if (length(dtheta) != 1)
       stop (ntheta,", should have length 1")
  }

  checkDim2 <- function (dr,dtheta,nx="dr",ntheta="dtheta") {
    if (length(dr) != dimC[1]+1)
      stop (nx,", should have length equal to dim(C)[1]+1 = ",dimC[1]+1," it is ",length(dr))
 
    if (length(dtheta) != dimC[2]+1)
       stop (ntheta,", should have length equal to dim(C)[2]+1 = ",dimC[2]+1," it is ",length(dtheta))
  }

  dimC <- dim (C) 
  Type <- 2
  
  if (max(theta)> 2*pi)
     stop ("theta should be < 2pi")
  if (min(theta)< 0)
     stop ("theta should be > 0 ")
  
  checkDim2(r, theta,"r, the grid in x-direction",
                   "theta, the grid in theta-direction") 
  # central values
  Nr       <- dimC[1]
  r.c      <- 0.5*(r[-1]+r[-(Nr+1)])
  dr       <- diff(r)
  drint    <- c(r.c[1]-r[1],diff(r.c),r[Nr+1]-r.c[Nr])
  divr     <- 1/r
  divr[is.na(divr)] <- 0
  
  Ntheta     <- dimC[2]
  theta.c    <- 0.5*(theta[-1]+theta[-(Ntheta +1)])
  dtheta     <- diff(theta)
  dthetaint  <- c(theta.c[1]-theta[1],diff(theta.c),theta[Ntheta+1]-theta.c[Ntheta])

  checkDim (D.r,D.theta,"D.r, the diffusion in x-direction ",
                      "D.theta, the diffusion in theta-direction ")
  
# check the dimensionality of the boundaries  
    if (! is.null(cyclicBnd)) {
# check if other boundaries not prescribed in this direction
    }
    if (! is.null(C.r.up)){ 
     if( length(C.r.up) != 1 && length(C.r.up) != Ntheta)
       stop ("'C.r.up' should have length 1 or equal to ", Ntheta)
    } else if (! is.null(flux.r.up)) { 
     if( length(flux.r.up) != 1 && length(flux.r.up) != Ntheta)
      stop ("'flux.r.up' should have length 1 or equal to ", Ntheta)
    } else if (!1 %in% cyclicBnd) 
      C.r.up = C[1,]

    if (! is.null(C.r.down)) {
     if( length(C.r.down) != 1 && length(C.r.down) != Ntheta)
      stop ("'C.r.down' should have length 1 or equal to ", Ntheta)
    } else if (! is.null(flux.r.down)){
     if(length(flux.r.down) != 1 && length(flux.r.down) != Ntheta)
      stop ("'flux.r.down' should have length 1 or equal to ", Ntheta)
    } else if (!1 %in% cyclicBnd)
      C.r.down = C[Nr,]

    if (! is.null(C.theta.up)) {
     if ( length(C.theta.up) != 1 && length(C.theta.up) != Nr)
      stop ("'C.theta.up' should have length 1 or equal to ", Nr)
    } else if (! is.null(flux.theta.up)){
     if (length(flux.theta.up) != 1 && length(flux.theta.up) != Nr)
       stop ("'flux.theta.up' should have length 1 or equal to ", Nr)
    } else if (!2 %in% cyclicBnd) 
      stop ("'flux.theta.up' OR 'C.theta.up' should be specified")

    if (! is.null(C.theta.down)) { 
     if(length(C.theta.down) != 1 && length(C.theta.down) != Nr)
      stop ("'C.theta.down' should have length 1 or equal to ", Nr)
    } else if (! is.null(flux.theta.down)) {
      if( length(flux.theta.down) != 1 && length(flux.theta.down) != Nr)
       stop ("'flux.theta.down' should have length 1 or equal to ", Nr)
    } else if (!2 %in% cyclicBnd) 
      stop ("'flux.theta.down' OR 'C.theta.down' should be specified")


# ------------------------------------------------------------------------------
# Function body: calculation
# ------------------------------------------------------------------------------

## Initialise 'fluxes' in all directions
  Flux.r <- 0
  Flux.theta <- 0

## First rewrite boundary values as "concentration"
## then perform diffusive transport
   if (1 %in% cyclicBnd) {
     C.r.up   <- (C[1,]*drint[Nr+1] + C[Nr,]*drint[1]) /(drint[1]+drint[Nr+1])
     C.r.down <- C.r.up
   }
   if (is.null(C.r.up))
     C.r.up <-  flux.r.up / D.r * drint[1]+ C[1,]
   if (is.null(C.r.down))
     C.r.down <- - flux.r.down / D.r * drint[Nr+1] + C[Nr,]
  
   Flux.r <- D.r * (rbind(C.r.up,C, deparse.level = 0) -
                    rbind(C,C.r.down, deparse.level = 0))/drint  

   if (2 %in% cyclicBnd) {
     C.theta.up   <- (C[,1]*dthetaint[Ntheta+1] + C[,Ntheta]*dthetaint[1]) /
        (dthetaint[1]+dthetaint[Ntheta+1])
     C.theta.down <- C.theta.up
   }
   if (is.null(C.theta.up))
     C.theta.up <- flux.theta.up/D.theta * dthetaint[1]*r.c + C[,1]
   if (is.null(C.theta.down))
     C.theta.down <- - flux.theta.down /D.theta * dthetaint[Ntheta+1]*r.c + C[,Ntheta]

   
   Flux.theta <- D.theta * (cbind(C.theta.up,C,deparse.level = 0) - 
                        cbind(C,C.theta.down, deparse.level = 0))/
                matrix(data= dthetaint,nrow=Nr,ncol=(Ntheta+1),byrow=TRUE) /r.c

## Calculate rate of change = flux gradient 

  dFdtheta <- 0.
  dFdr   <- -diff(Flux.r * r)/dr/r.c
  dFdtheta <- -t(diff(t(Flux.theta))/dtheta)/r.c

 if (full.output)
  return (list (dC = dFdr + dFdtheta ,                  # Rate of change
                C.r.up     = C.r.up,
                C.r.down   = C.r.down,
                C.theta.up   = C.theta.up,
                C.theta.down = C.theta.down,
                r.flux   = Flux.r,                 
                theta.flux = Flux.theta,
                flux.r.up = Flux.r[1,],
                flux.r.down = Flux.r[Nr+1,],
                flux.theta.up = Flux.theta[,1],
                flux.theta.down = Flux.theta[,Ntheta+1]
                
                ))  
 else
  return (list (dC = dFdr + dFdtheta ,                 
                flux.r.up = Flux.r[1,],
                flux.r.down = Flux.r[Nr+1,],
                flux.theta.up = Flux.theta[,1],
                flux.theta.down = Flux.theta[,Ntheta+1]
                ))  
 
} # end tran.polar
