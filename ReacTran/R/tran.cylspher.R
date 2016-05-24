##==============================================================================
## Diffusive transport in cylindrical coordinates (r, theta, z)
##==============================================================================

tran.cylindrical <- function(C, C.r.up = NULL, C.r.down = NULL, 
                             C.theta.up = NULL, C.theta.down = NULL, 
                             C.z.up = NULL, C.z.down = NULL, 
                             flux.r.up = NULL, flux.r.down = NULL, 
                             flux.theta.up = NULL, flux.theta.down = NULL,          
                             flux.z.up = NULL, flux.z.down = NULL, 
                             cyclicBnd = NULL,
                             D.r = NULL, D.theta = D.r, D.z = D.r, 
                             r = NULL,  theta = NULL, z = NULL) 
{
# ------------------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------------------

  dimens <- dim(C)
  if (length(dimens) != 3)
    stop("'C' should an array of dimensionality 3")
  
  Nr     <- dimens[1]
  Ntheta <- dimens[2]
  Nz     <- dimens[3]

  if (length(r) != Nr+1)
    stop("Length of 'r' should equal first dimension of 'C' + 1")
  if (length(theta) != Ntheta+1)
    stop("Length of 'theta' should equal second dimension of 'C' + 1")
  if (length(z) != Nz+1)
    stop("Length of 'z' should equal third dimension of 'C' + 1")

  if (max(theta) > 2 * pi)                                                               
    stop("theta should be < 2pi")                                                      
  if (min(theta) < 0)                                                                    
    stop("theta should be > 0 ")                                                       

# boundary conditions default = zero-gradient
  Bc.r.Up       <- 3
  Bc.r.Down     <- 3
  Bc.theta.Up   <- 3
  Bc.theta.Down <- 3
  Bc.z.Up       <- 3
  Bc.z.Down     <- 3

# boundaries in r-direction
  if (! is.null(C.r.up))
    Bc.r.Up <- 2
  else
    C.r.up <- 0.  

  if (! is.null(flux.r.up))
    Bc.r.Up <- 1
  else
    flux.r.up <- 0.  
    
  if (! is.null(C.r.down))
    Bc.r.Down <- 2
  else
    C.r.down <- 0.  

  if (! is.null(flux.r.down))
    Bc.r.Down <- 1
  else            
    flux.r.down <- 0.

  if (1 %in% cyclicBnd) {
    Bc.r.Up   <- 5
    Bc.r.Down <- 5
  }  

# boundaries in theta-direction
  if (! is.null(C.theta.up))
    Bc.theta.Up <- 2
  else
    C.theta.up <- 0.  

  if (! is.null(flux.theta.up))
    Bc.theta.Up <- 1
  else
    flux.theta.up <- 0.  
    
  if (! is.null(C.theta.down))
    Bc.theta.Down <- 2
  else
    C.theta.down <- 0.  

  if (! is.null(flux.theta.down))
    Bc.theta.Down <- 1
  else            
    flux.theta.down <- 0.

  if (2 %in% cyclicBnd) {
    Bc.theta.Up   <- 5
    Bc.theta.Down <- 5
  }  

# boundaries in z-direction
  if (! is.null(C.z.up))
    Bc.z.Up <- 2
  else
    C.z.up <- 0.  

  if (! is.null(flux.z.up))
    Bc.z.Up <- 1
  else
    flux.z.up <- 0.  
    
  if (! is.null(C.z.down))
    Bc.z.Down <- 2
  else
    C.z.down <- 0.  

  if (! is.null(flux.z.down))
    Bc.z.Down <- 1
  else            
    flux.z.down <- 0.

  if (3 %in% cyclicBnd) {
    Bc.phi.Up   <- 5
    Bc.phi.Down <- 5
  }  

# Diffusion coefficients
  if (! is.null(D.r))
    D.r <- rep(D.r, length.out=Nr+1)
  if (! is.null(D.theta))
    D.theta <- rep(D.theta, length.out=Ntheta+1)
  if (! is.null(D.z))
    D.z <- rep(D.z, length.out=Nz+1)

  if (length(C.r.up) + length(flux.r.up) + 
      length(C.r.down) + length(flux.r.down) +
      length(C.theta.up) + length(flux.theta.up)+
      length(C.theta.down) + length(flux.theta.down)+
      length(C.z.up) + length(flux.z.up) + 
      length(C.z.down)+length(flux.z.down)!=12)
   stop("length of all 'boundary conditions' should be 1") 


    tr <- .Fortran("diffcyl", as.integer(Nr), as.integer(Ntheta), 
      as.integer(Nz), as.double(C),  
      as.double(C.r.up), as.double(C.r.down),  
      as.double(C.theta.up), as.double(C.theta.down),
      as.double(C.z.up), as.double(C.z.down),
      as.double(flux.r.up),  as.double(flux.r.down), 
      as.double(flux.theta.up),  as.double(flux.theta.down), 
      as.double(flux.z.up),  as.double(flux.z.down), 
      as.integer(Bc.r.Up), as.integer(Bc.r.Down),  
      as.integer(Bc.theta.Up), as.integer(Bc.theta.Down),  
      as.integer(Bc.z.Up), as.integer(Bc.z.Down),  
      as.double(D.r), as.double(D.theta), as.double(D.z), 
      as.double(r), as.double(theta), as.double(z),
      Fluxrup= numeric(Ntheta*Nz),Fluxrdown= numeric(Ntheta*Nz),
      Fluxthetaup= numeric(Nr*Nz),Fluxthetadown= numeric(Nr*Nz),
      FluxZup= numeric(Nr*Ntheta),FluxZdown= numeric(Nr*Ntheta),
      dC = numeric(Nr*Ntheta*Nz))  

  list(dC = array(dim=c(Nr, Ntheta, Nz), tr$dC),
       flux.r.up       = matrix(nrow=Ntheta, ncol = Nz, tr$Fluxrup),
       flux.r.down     = matrix(nrow=Ntheta, ncol = Nz, tr$Fluxrdown),
       flux.theta.up   = matrix(nrow=Nr, ncol = Nz, tr$Fluxthetaup),
       flux.theta.down = matrix(nrow=Nr, ncol = Nz, tr$Fluxthetadown),
       flux.z.up       = matrix(nrow=Nr, ncol = Ntheta, tr$FluxZup),
       flux.z.down     = matrix(nrow=Nr, ncol = Ntheta, tr$FluxZdown)
       )
}   


##==============================================================================
## Diffusive transport in spherical coordinates (r, theta, phi)
##==============================================================================

tran.spherical <- function(C, C.r.up = NULL, C.r.down = NULL, 
                             C.theta.up = NULL, C.theta.down = NULL, 
                             C.phi.up = NULL, C.phi.down = NULL, 
                             flux.r.up = NULL, flux.r.down = NULL, 
                             flux.theta.up = NULL, flux.theta.down = NULL,          
                             flux.phi.up = NULL, flux.phi.down = NULL, 
                             cyclicBnd = NULL,
                             D.r = NULL, D.theta = D.r, D.phi = D.r, 
                             r = NULL,  theta = NULL, phi = NULL) 
{
# ------------------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------------------

  dimens <- dim(C)
  if (length(dimens) != 3)
    stop("'C' should an array of dimensionality 3")
  
  Nr     <- dimens[1]
  Ntheta <- dimens[2]
  Nphi     <- dimens[3]

  if (length(r) != Nr+1)
    stop("Length of 'r' should equal first dimension of 'C' + 1")
  if (length(theta) != Ntheta+1)
    stop("Length of 'theta' should equal second dimension of 'C' + 1")
  if (length(phi) != Nphi+1)
    stop("Length of 'phi' should equal third dimension of 'C' + 1")

  if (max(theta) > 2 * pi)                                                               
    stop("theta should be <= 2*pi")                                                      
  if (min(theta) < 0)                                                                    
    stop("theta should be >= 0 ")                                                       

  if (max(phi) > 2 * pi)                                                               
    stop("phi should be <= 2*pi")                                                      
  if (min(phi) < 0)                                                                    
    stop("phi should be >= 0 ")                                                       

# boundary conditions default = zero-gradient
  Bc.r.Up       <- 3
  Bc.r.Down     <- 3
  Bc.theta.Up   <- 3
  Bc.theta.Down <- 3
  Bc.phi.Up       <- 3
  Bc.phi.Down     <- 3

# boundaries in r-direction
  if (! is.null(C.r.up))
    Bc.r.Up <- 2
  else
    C.r.up <- 0.  

  if (! is.null(flux.r.up))
    Bc.r.Up <- 1
  else
    flux.r.up <- 0.  
    
  if (! is.null(C.r.down))
    Bc.r.Down <- 2
  else
    C.r.down <- 0.  

  if (! is.null(flux.r.down))
    Bc.r.Down <- 1
  else            
    flux.r.down <- 0.

  if (1 %in% cyclicBnd) {
    Bc.r.Up   <- 5
    Bc.r.Down <- 5
  }  

# boundaries in theta-direction
  if (! is.null(C.theta.up))
    Bc.theta.Up <- 2
  else
    C.theta.up <- 0.  

  if (! is.null(flux.theta.up))
    Bc.theta.Up <- 1
  else
    flux.theta.up <- 0.  
    
  if (! is.null(C.theta.down))
    Bc.theta.Down <- 2
  else
    C.theta.down <- 0.  

  if (! is.null(flux.theta.down))
    Bc.theta.Down <- 1
  else            
    flux.theta.down <- 0.

  if (2 %in% cyclicBnd) {
    Bc.theta.Up   <- 5
    Bc.theta.Down <- 5
  }  

# boundaries in phi-direction
  if (! is.null(C.phi.up))
    Bc.phi.Up <- 2
  else
    C.phi.up <- 0.  

  if (! is.null(flux.phi.up))
    Bc.phi.Up <- 1
  else
    flux.phi.up <- 0.  
    
  if (! is.null(C.phi.down))
    Bc.phi.Down <- 2
  else
    C.phi.down <- 0.  

  if (! is.null(flux.phi.down))
    Bc.phi.Down <- 1
  else            
    flux.phi.down <- 0.

  if (3 %in% cyclicBnd) {
    Bc.phi.Up   <- 5
    Bc.phi.Down <- 5
  }  

# Diffusion coefficients
  if (! is.null(D.r))
    D.r <- rep(D.r, length.out=Nr+1)
  if (! is.null(D.theta))
    D.theta <- rep(D.theta, length.out=Ntheta+1)
  if (! is.null(D.phi))
    D.phi <- rep(D.phi, length.out=Nphi+1)

  if (length(C.r.up) + length(flux.r.up) + 
      length(C.r.down) + length(flux.r.down) +
      length(C.theta.up) + length(flux.theta.up)+
      length(C.theta.down) + length(flux.theta.down)+
      length(C.phi.up) + length(flux.phi.up) + 
      length(C.phi.down)+length(flux.phi.down)!=12)
   stop("length of all 'boundary conditions' should be 1") 


    tr <- .Fortran("diffspher", as.integer(Nr), as.integer(Ntheta), 
      as.integer(Nphi), as.double(C),  
      as.double(C.r.up), as.double(C.r.down),  
      as.double(C.theta.up), as.double(C.theta.down),
      as.double(C.phi.up), as.double(C.phi.down),
      as.double(flux.r.up),  as.double(flux.r.down), 
      as.double(flux.theta.up),  as.double(flux.theta.down), 
      as.double(flux.phi.up),  as.double(flux.phi.down), 
      as.integer(Bc.r.Up), as.integer(Bc.r.Down),  
      as.integer(Bc.theta.Up), as.integer(Bc.theta.Down),  
      as.integer(Bc.phi.Up), as.integer(Bc.phi.Down),  
      as.double(D.r), as.double(D.theta), as.double(D.phi), 
      as.double(r), as.double(theta), as.double(phi),
      Fluxrup= numeric(Ntheta*Nphi),Fluxrdown= numeric(Ntheta*Nphi),
      Fluxthetaup= numeric(Nr*Nphi),Fluxthetadown= numeric(Nr*Nphi),
      Fluxphiup= numeric(Nr*Ntheta),Fluxphidown= numeric(Nr*Ntheta),
      dC = numeric(Nr*Ntheta*Nphi))  

  list(dC = array(dim=c(Nr, Ntheta, Nphi), tr$dC),
       flux.r.up       = matrix(nrow=Ntheta, ncol = Nphi, tr$Fluxrup),
       flux.r.down     = matrix(nrow=Ntheta, ncol = Nphi, tr$Fluxrdown),
       flux.theta.up   = matrix(nrow=Nr, ncol = Nphi, tr$Fluxthetaup),
       flux.theta.down = matrix(nrow=Nr, ncol = Nphi, tr$Fluxthetadown),
       flux.phi.up     = matrix(nrow=Nr, ncol = Ntheta, tr$Fluxphiup),
       flux.phi.down   = matrix(nrow=Nr, ncol = Ntheta, tr$Fluxphidown)
       )
}   
