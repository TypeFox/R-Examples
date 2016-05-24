## Calculates Mrad and Mtot for a dipole near a multilayer, using the angular decomposition of the dipole field into plane waves
## Mrad is done in one step, but the integration for Mtot is divided in 3 regions
## thus we define 4 different integrands


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the radiative dipole decay rates near a multilayer interface. 
##' @title integrand_rad
##' @export
##' @param d distance in nm
##' @param angle angle in radians
##' @param wavelength wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param GL logical: result formatted for use with Gauss Legendre quadrature
##' @author baptiste Auguie
##' @family integrands dipole
integrand_rad <- function(d = 10, angle, wavelength,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0),  GL=FALSE){

  ## for 0 < q < 1, i.e. 0 < u < 1
  
  ## define constants
  k0 <- 2*pi/wavelength
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Ntheta <- length(angle)
  
  cost <- cos(angle)
  sint <- sin(angle)
  
  rp <- -1*recursive_fresnelcpp(wavelength=wavelength,
                           q = sint,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(wavelength=wavelength,
                           q = sint,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection

  phase <- exp(2i*d*outer(k1,cost))
  
  integrand.p <- Mod(matrix(1, Nlambda, Ntheta, byrow=TRUE) + rp * phase)^2 *
    matrix(sint^3, Nlambda, Ntheta, byrow=TRUE)
  
  integrand.s <- (Mod(matrix(1, Nlambda, Ntheta, byrow=TRUE) + rs * phase)^2 +
                  Mod(matrix(1, Nlambda, Ntheta, byrow=TRUE) - rp * phase)^2 *
                  matrix(cost^2, Nlambda, Ntheta, byrow=TRUE)) *
                    matrix(sint, Nlambda, Ntheta, byrow=TRUE)

  if(GL)
    return(list(integrand.p = integrand.p, integrand.s = integrand.s)) else 
  rbind(integrand.p, integrand.s)
}


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface. Transformed part I1 (radiative)
##' from u=0 to 1
##' @title integrand_nr1
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(1-q^2)
##' @param wavelength wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param GL logical: result formatted for use with Gauss Legendre quadrature
##' @author baptiste Auguie
##' @family integrands dipole
integrand_nr1 <- function(d=10, u, wavelength,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0), GL=FALSE){

  ## integrand1 is for 0 < q < 1, i.e. 0 < u < 1
  
  ## define constants
  k0 <- 2*pi/wavelength
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  rp <- -1*recursive_fresnelcpp(wavelength=wavelength,
                           q = sqrt(1 - u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(wavelength=wavelength,
                           q = sqrt(1 - u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
    phase <- exp(2i*d*outer(k1,u))
  
    integrand.p <- Re(matrix(1 - u^2, Nlambda, Nq, byrow=TRUE) * rp * phase)
    integrand.s <- Re(( rs - rp*matrix(u^2, Nlambda, Nq, byrow=TRUE)) * phase)
      
  if(GL)
    return(list(integrand.p = integrand.p, integrand.s = integrand.s)) else 
  rbind(integrand.p, integrand.s)
}

##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface. Transformed part I2
##' from u=0 to ucut
##' @title integrand_nr2
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(q^2 - 1)
##' @param wavelength wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param GL logical: result formatted for use with Gauss Legendre quadrature
##' @author baptiste Auguie
##' @family integrands dipole
integrand_nr2 <- function(d=10, u, wavelength,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0),  GL=FALSE){

  ## integrand2 is for 1 < q < infty, i.e. 0 < u < infty

  ## define constants
  k0 <- 2*pi/wavelength
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  rp <- -1*recursive_fresnelcpp(wavelength=wavelength,
                           q = sqrt(1 + u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(wavelength=wavelength,
                           q = sqrt(1 + u^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
  ## phase is now purely real
  phase <- exp(-2*d*outer(k1,u))

  ## take the imaginary part now
  integrand.p <- matrix(1 + u^2, Nlambda, Nq, byrow=TRUE) * Im( rp ) * phase
  integrand.s <- Im(rs + rp*matrix(u^2, Nlambda, Nq, byrow=TRUE)) * phase
  
  if(GL)
    return(list(integrand.p = integrand.p, integrand.s = integrand.s)) else 
  rbind(integrand.p, integrand.s)
}


##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface. Transformed part III
##' from u=ucut to infinity
##' @title integrand_nr3
##' @export
##' @param d distance in nm
##' @param u transformed normalised in-plane wavevector sqrt(q^2 - 1)
##' @param ucut limit of the integral
##' @param wavelength wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param GL logical: result formatted for use with Gauss Legendre quadrature
##' @author baptiste Auguie
##' @family integrands dipole
integrand_nr3 <- function(d=10, u, ucut, wavelength,
                       epsilon = list(incident=1.5^2, 1.0^2),
                       thickness = c(0, 0), GL=FALSE){

  ## define constants
  k0 <- 2*pi/wavelength
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(u)
  
  ## integrand2 is for ucut < u < infty
  ## performing a change of variable mapping u in [ucut, infty) -> [0,1]
  ## change of variables
  ## \int_a^\infty f(x)dx = \int_0^1 f(a + t/(1-t)). 1 / (1-t)^2 dt
  ## as suggested on http://ab-initio.mit.edu/wiki/index.php/Cubature

  ## new variable
  t <-  ucut + u / (1 - u)
  ## Jacobian of transformation
  Jac <-  matrix(1 / (1 - u)^2, Nlambda, Nq, byrow=TRUE)
  
  
  rp <- -1*recursive_fresnelcpp(wavelength=wavelength,
                           q = sqrt(1 + t^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(wavelength=wavelength,
                           q = sqrt(1 + t^2),
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
  ## phase is now purely real
  phase <- exp(-2*d*outer(k1,t))

  ## take the imaginary part now
  integrand.p <- matrix(1 + t^2, Nlambda, Nq, byrow=TRUE) * Im( rp ) * phase * Jac
  integrand.s <- Im(rs + rp*matrix(t^2, Nlambda, Nq, byrow=TRUE)) * phase * Jac
  
  if(GL)
    return(list(integrand.p = integrand.p, integrand.s = integrand.s)) else 
  rbind(integrand.p, integrand.s)
}

##' Dipole decay rates near a multilayer interface
##'
##' dipole decay rates near a multilayer interface
##' @title dipole
##' @export
##' @param d distance in nm
##' @param wavelength wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param qcut transition between regions 2 and 3
##' @param rel.err relative error
##' @param Nquadrature1 maximum number of quadrature points in radiative region
##' @param Nquadrature2 maximum number of quadrature points in SPPs region
##' @param Nquadrature3 maximum number of quadrature points in dipole image region
##' @param GL logical: use Gauss Legendre quadrature,  or cubature::adaptIntegrate
##' @param show.messages logical, display integration info
##' @author baptiste Auguie
dipole <- function(d=1,
                   wavelength,
                   epsilon = list(incident=1.0^2),
                   thickness = c(0, 0), qcut=NULL, rel.err = 1e-3,
                   Nquadrature1 = 1e3, Nquadrature2 = 1e4, Nquadrature3 = 1e4,
                   GL = FALSE, 
                   show.messages=TRUE){
   
  Nlambda <- length(wavelength)
  

  if(GL){
    GL1 <- gauss.quad(Nquadrature1)
    GL2 <- gauss.quad(Nquadrature2)
    GL3 <- gauss.quad(Nquadrature3)
  } 

  ## if no qcut provided, estimate one from max of
  ## all possible SPP dispersions
  if(is.null(qcut)){
    qcut <- 1.1

    epsilon_norm <- do.call(cbind, epsilon)
    
    for(ii in seq(1, length(epsilon) - 1)){
      qspp <- sqrt(epsilon_norm[,ii] / epsilon_norm[,1])*
        sqrt(epsilon_norm[,ii+1] / (epsilon_norm[,ii] + epsilon_norm[,ii+1]))
      
      qcut <- max(qcut, max(Re(qspp)))
    }

    if(show.messages)
      message(paste("using qcut=", round(qcut,2)))
    
  }
  
  ## integration from 0 to 1 for the transformed radiative bit

  if(GL){
    
    umax1 <- 1; umin1 <- 0;
    C1 <- (umax1 - umin1)/2 ; D1 <- (umax1+umin1)/2
    unodes1 <- C1 * GL1$nodes + D1
    uweights1 <- GL1$weights * C1
    
    Nu1 <- length(unodes1)
    
#     message("here")
    in1 <- integrand_nr1(u=unodes1,
                        d=d, wavelength=wavelength,
                        epsilon=epsilon, thickness=thickness, GL=TRUE)
    weights1 <- matrix(uweights1, nrow=Nlambda, ncol=Nu1, byrow=TRUE)
    
    integral1.perp <- rowSums(in1$integrand.p*weights1)
    integral1.par <- rowSums(in1$integrand.s*weights1)
    
  } else {    
    
    in1 <- adaptIntegrate(integrand_nr1, lowerLimit = 0,
                          upperLimit = 1,
                          fDim = 2*Nlambda, tol=rel.err,
                          maxEval = Nquadrature1, 
                          d=d, wavelength=wavelength,
                          epsilon=epsilon, thickness=thickness)
    
    integral1.perp <- in1$integral[seq(1,Nlambda)]
    integral1.par <- in1$integral[seq(Nlambda+1,2*Nlambda)]
    
  }
  
  ## integration from 0 to ucut
  
  if(GL){
    
    ucut <- sqrt(qcut^2 - 1)
    
    umax2 <- ucut; umin2 <- 0;
    C2 <- (umax2 - umin2)/2 ; D2 <- (umax2 + umin2)/2
    
    unodes2 <- C2 * GL2$nodes + D2
    uweights2 <- GL2$weights * C2 
    
    Nu2 <- length(unodes2)
    
    in2 <- integrand_nr2(u=unodes2,
                        d=d, wavelength=wavelength,
                    epsilon=epsilon, thickness=thickness,  GL=TRUE)
    
    weights2 <- matrix(uweights2, nrow=Nlambda, ncol=Nu2, byrow=TRUE)
    
    integral2.perp <- rowSums(in2$integrand.p*weights2)
    integral2.par <- rowSums(in2$integrand.s*weights2)
     
  } else {
    
    ucut <- sqrt(qcut^2 - 1)
    
    in2 <- adaptIntegrate(integrand_nr2, lowerLimit = 0,
                          upperLimit = ucut,
                          fDim = 2*Nlambda, tol=rel.err,
                          maxEval = Nquadrature2,
                          d=d, wavelength=wavelength,
                          epsilon=epsilon, thickness=thickness)
    
    integral2.perp <- in2$integral[seq(1,Nlambda)]
    integral2.par <- in2$integral[seq(Nlambda+1,2*Nlambda)]
  }
  
  ## integration from ucut to Inf
  ## integrand performing a change of variable mapping u in [ucut, infty) -> t in [0,1]
  
  if(GL){
    ## performing a change of variable mapping u in [ucut, infty) -> [0,1]
    ## change of variables
    ## \int_a^\infty f(x)dx = \int_0^1 f(a + t/(1-t)). 1 / (1-t)^2 dt
    ## as suggested on http://ab-initio.mit.edu/wiki/index.php/Cubature

    umax3 <- 1; umin3 <- 0;
    C3 <- (umax3 - umin3)/2 ; D3 <- (umax3 + umin3)/2
    
    unodes3 <- C3 * GL3$nodes + D3
    uweights3 <- GL3$weights * C3 * 1 / (1 - unodes3)^2
    unodes3 <- ucut + unodes3 / (1 - unodes3)
    
    Nu3 <- length(unodes3)
    
    in3 <- integrand_nr2(u=unodes3,
                        d=d, wavelength=wavelength,
                        epsilon=epsilon, thickness=thickness,  GL=TRUE)
    
    weights3 <- matrix(uweights3, nrow=Nlambda, ncol=Nu3, byrow=TRUE)
    
    integral3.perp <- rowSums(in3$integrand.p*weights3)
    integral3.par <- rowSums(in3$integrand.s*weights3)
    
  } else {
    # use integrand_nr2 as we do the variable transformation outside (? to check)
    in3 <- adaptIntegrate(integrand_nr3, lowerLimit = 0,
                          upperLimit = 1,
                          fDim = 2*Nlambda, tol=rel.err,
                          maxEval = Nquadrature3, 
                          ucut=ucut, d=d, wavelength=wavelength,
                          epsilon=epsilon, thickness=thickness)
    
    integral3.perp <- in3$integral[seq(1,Nlambda)]
    integral3.par <- in3$integral[seq(Nlambda+1,2*Nlambda)]
    
  }
  ## Mrad
  
  if(GL){
    ## for Mrad, we use the same integration points as GL1 because we study the radiative region
    
    anglemax <- pi/2; anglemin <- 0;
    C4 <- (anglemax - anglemin)/2 ; D4 <- (anglemax + anglemin)/2
    anglenodes <- C4 * GL1$nodes + D4
    angleweights <- GL1$weights * C4
    
    Ntheta <- length(anglenodes)
    
    in4 <- integrand_rad(angle=anglenodes,
                           d=d, wavelength=wavelength,
                           epsilon=epsilon, thickness=thickness, GL=TRUE)
    
    weights4 <- matrix(angleweights, nrow=Nlambda, ncol=Ntheta, byrow=TRUE)
    
    Mrad.perp <- 3/4 * rowSums(in4$integrand.p*weights4)
    Mrad.par <- 3/8 * rowSums(in4$integrand.s*weights4)
    
  } else {
    
    in4 <- adaptIntegrate(integrand_rad, lowerLimit = 0,
                          upperLimit = pi/2,
                          fDim = 2*Nlambda, tol=rel.err,
                          maxEval = Nquadrature1, 
                          d=d, wavelength=wavelength,
                          epsilon=epsilon, thickness=thickness)
    
    Mrad.perp <-  3/4 *in4$integral[seq(1,Nlambda)]
    Mrad.par <-  3/8 *in4$integral[seq(Nlambda+1,2*Nlambda)]
    
    evaluations <- c(in1$functionEvaluations,
                     in2$functionEvaluations,
                     in3$functionEvaluations,
                     in4$functionEvaluations)
    
    integration <- sprintf("relative integration errors were: %.3e for I1,  %.3e for I2,  %.3e for I3,  %.3e for I4; with %i, %i, %i, %i respective function evaluations.\n",
            max(in1$error), max(in2$error), max(in3$error), max(in4$error), 
            evaluations[1],evaluations[2],evaluations[3],evaluations[4])
  }
  
  
  results <- data.frame(wavelength=wavelength,
                  Mtot.perp = 1 + 3/2*(integral1.perp + integral2.perp + integral3.perp),
                  Mtot.par = 1 + 3/4*(integral1.par + integral2.par + integral3.par),
                  Mrad.perp =  Mrad.perp, Mrad.par = Mrad.par)
  if(!GL && show.messages)
    message(integration)
  if(!GL)
  comment(results) <- integration

  invisible(results)
}

