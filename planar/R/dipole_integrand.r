
##' Total decay rate of a dipole near a multilayer interface
##'
##' Integrand without transformation of variables
##' @title integrand_mtot
##' @export
##' @param d distance in nm
##' @param q normalised in-plane wavevector in [0, infty) 
##' @param wavelength wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @family integrands dipole
##' @author baptiste Auguie
integrand_mtot <- function(d=10, q, wavelength,
                             epsilon = list(incident=1.5^2, 1.0^2),
                             thickness = c(0, 0)){
  
  ## define constants
  k0 <- 2*pi/wavelength
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(q)
  
  u <- sqrt(1 - q^2 + 0i)
  rp <- -1* recursive_fresnelcpp(wavelength=wavelength,
                                 q = q,
                                 epsilon=epsilon,
                                 thickness=thickness,
                                 polarisation="p")$reflection
  
  rs <- recursive_fresnelcpp(wavelength=wavelength,
                           q = q,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
  phase <- exp(2i*d*outer(k1,u))
  
  integrand.p <- Re(matrix(q^3 / u, Nlambda, Nq, byrow=TRUE) * rp * phase)
  integrand.s <- Re( (rs / matrix(u, Nlambda, Nq, byrow=TRUE) -
                      rp * matrix(u, Nlambda, Nq, byrow=TRUE)) *
                    matrix(q, Nlambda, Nq, byrow=TRUE) * phase)
   
  list(integrand.p = integrand.p, integrand.s = integrand.s)
}


##' Dipole total decay rate near a multilayer interface
##'
##' direct application of the textbook formula using integrand_mtot; performs poorly compared to the transformed version in \code{dipole}
##' @title dipole_direct
##' @export
##' @param d distance in nm
##' @param wavelength wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param Nquadrature1 quadrature points in radiative region
##' @param Nquadrature2 quadrature points in SPPs region
##' @param Nquadrature3 quadrature points in dipole image region
##' @param qcut transition between regions 2 and 3
##' @param qmax maximum q of region 3
##' @param show.messages logical, display integration info
##' @family dipole
##' @author baptiste Auguie
dipole_direct <- function(d=1,
                   wavelength ,
                   epsilon = list(incident=1.0^2),
                   thickness = c(0, 0),
                   Nquadrature1 = 50, Nquadrature2 = 200, Nquadrature3 = 50,
                   qcut = NULL, qmax = Inf, show.messages=TRUE){

  GL1 <- gauss.quad(Nquadrature1)
  GL2 <- gauss.quad(Nquadrature2)
  GL3 <- gauss.quad(Nquadrature3)

  Nq1 <- length(GL1$nodes)
  Nq2 <- length(GL2$nodes)
  Nq3 <- length(GL3$nodes)
  
  Nlambda <- length(wavelength)
  
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
      message(paste("using qcut=",round(qcut,2)))
    
  }

  ## integration from 0 to 1
  qmax1 <- 1; qmin1 <- 0;
  C1 <- (qmax1 - qmin1)/2 ; D1 <- (qmax1+qmin1)/2
  qnodes1 <- C1 * GL1$nodes + D1
  qweights1 <- GL1$weights * C1
  
  in1 <- integrand_mtot(q=qnodes1,
                          d=d, wavelength=wavelength,
                          epsilon=epsilon, thickness=thickness)
      
  weights1 <- matrix(qweights1, nrow=Nlambda, ncol=Nq1, byrow=TRUE)

  integral1.perp <- rowSums(in1$integrand.p*weights1)
  integral1.par <- rowSums(in1$integrand.s*weights1)
  
  
  ## integration from 1 to qcut
  qmax2 <- qcut; qmin2 <- 1;
  C2 <- (qmax2 - qmin2)/2 ; D2 <- (qmax2+qmin2)/2
  qnodes2 <- C2 * GL2$nodes + D2
  qweights2 <- GL2$weights * C2
  
  in2 <- integrand_mtot(q=qnodes2,
                                  d=d, wavelength=wavelength,
                                  epsilon=epsilon, thickness=thickness)
  
  weights2 <- matrix(qweights2, nrow=Nlambda, ncol=Nq2, byrow=TRUE)
  
  integral2.perp <- rowSums(in2$integrand.p*weights2)
  integral2.par <- rowSums(in2$integrand.s*weights2)
  

  ## integration from qcut to qmax
  if(is.finite(qmax)){
    ## straight integration from qcut to qmax
    qmax3 <- qmax; qmin3 <- qcut;
    C3 <- (qmax3 - qmin3)/2 ; D3 <- (qmax3+qmin3)/2
    
    qnodes3 <- C3 * GL3$nodes + D3
    qweights3 <- GL3$weights * C3
  
  } else {
    if(show.messages)
      message("performing a change of variable mapping [qcut, infty) -> [0,1]")
    ## change of variables
    ## \int_a^\infty f(x)dx = \int_0^1 f(a + t/(1-t)). 1 / (1-t)^2 dt
    ## as suggested on http://ab-initio.mit.edu/wiki/index.php/Cubature
    qmax3 <- 1; qmin3 <- 0;
    C3 <- (qmax3 - qmin3)/2 ; D3 <- (qmax3+qmin3)/2
    
    qnodes3 <- C3 * GL3$nodes + D3
    qweights3 <- GL3$weights * C3 * 1 / (1 - qnodes3)^2
    qnodes3 <- qcut + qnodes3 / (1 - qnodes3)
    
  }
  
  in3 <- integrand_mtot(q=qnodes3,
                          d=d, wavelength=wavelength,
                          epsilon=epsilon, thickness=thickness)
      
  weights3 <- matrix(qweights3, nrow=Nlambda, ncol=Nq3, byrow=TRUE)
  
  integral3.perp <- rowSums(in3$integrand.p*weights3)
  integral3.par <- rowSums(in3$integrand.s*weights3)
  
  data.frame(wavelength=wavelength,
             Mtot.perp = 1 + 3/2*(integral1.perp + integral2.perp + integral3.perp),
             Mtot.par = 1 + 3/4*(integral1.par + integral2.par + integral3.par) )
  
}
