## 
## Functions for circular dichroism (CD)
## 

integration_points <- function(method=c("cheap", "QMC", "GL", "grid"), 
                               Nquad = 30, init=TRUE){
  
  method <- match.arg(method)
  
  if(method == "cheap"){
    
    angles <- rbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                    c(pi/2, pi/2, 0), # +y is phi=pi/2, psi=0
                    c(pi/2, pi/2, pi/2)) # +z is phi=pi/2, psi=pi/2
    weights <- rep(1/nrow(angles), nrow(angles)) 
    return(list(angles=angles, weights=weights))
  }
  
  if(method == "QMC"){
    
    nodes <- randtoolbox::halton(Nquad, dim = 2, normal=FALSE, init=init)
    
    phi <- nodes[,1]*2*pi
    psi <- asin(2*nodes[,2] - 1)
    grid <- data.frame(phi=phi, theta=pi/2, psi=psi)
    weights <- rep(1/nrow(grid), nrow(grid))
    return(list(angles=grid, weights=weights))
  }
  
  if(method == "GL"){
    # scale the coordinates from (-1, 1) to (0, 2pi) and (-pi/2, pi/2) resp.
    phi1=0; phi2=2*pi; psi1=-pi/2; psi2=pi/2; 
    C1 = (phi2 - phi1) / 2;  D1 = (phi2 + phi1) / 2;
    C2 = (psi2 - psi1) / 2; D2 = (psi2 + psi1) / 2;
    
    rndN <- ceiling(sqrt(Nquad/2))
    GL_phi <- statmod::gauss.quad(2*rndN)
    GL_psi <- statmod::gauss.quad(rndN)
    
    phi = GL_phi$nodes*C1 + D1  
    psi = GL_psi$nodes*C2 + D2
    
    # grid of angles, theta is constant here
    grid <- expand.grid(phi=phi, theta=pi/2, psi=psi)
    # corresponding weights for 1D quadrature
    weights <- expand.grid(phi=GL_phi$weights, psi=GL_psi$weights)
    # combine the weigths for each point; cos(psi) comes from the Jacobian in the integral
    weights <- C1 * C2 / (4*pi) * cos(grid$psi) * weights$phi * weights$psi
    
    return(list(angles=grid, weights=weights))
  }  
  
  if(method == "grid"){
    ## remove end points that cause problems
    a <- seq(0 + 0.0001/Nquad,1-0.0001/Nquad, length=round(sqrt(Nquad))) 
    phi <- a*2*pi
    psi <- asin(2*a - 1)
    grid <- expand.grid(phi=phi, theta=pi/2, psi=psi)
    weights <- rep(1/nrow(grid), nrow(grid))
    return(list(angles=grid, weights=weights))
  }
  
}

##' Simulate a CD spectrum
##'
##' CD spectrum
##' @title circular_dichroism_spectrum 
##' @param cluster cluster (list)
##' @param material material
##' @param medium refractive index medium
##' @param Nquad number of integration points
##' @param averaging averaging method, using either Gauss Legendre quadrature (default), Quasi Monte Carlo, regular grid, or "cheap" (3 axes)
##' @param iterative logical, increase N until convergence (QMC only)
##' @param precision relative diff between two runs (QMC only)
##' @param Nmax maximum N if convergence not attained (QMC only)
##' @param dN iterative increase in N (QMC only)
##' @param full logical use full (retarded) dipolar field
##' @param progress print progress lines
##' @param verbose display messages
##' @param result.matrix logical return the results as a matrix
##' @importFrom randtoolbox halton
##' @importFrom statmod gauss.quad
##' @importFrom reshape2 melt
##' @export
##' @family user_level circular_dichroism
##' @author baptiste Auguie
##' @references
##' Y. Okada, Efficient numerical orientation averaging of light scattering properties with a quasi-Monte-Carlo method, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 109, Issue 9, June 2008, Pages 1719-1742.
circular_dichroism_spectrum <- function(cluster, material, medium=1.33, Nquad=100, 
                                        averaging = c("GL","QMC","grid", "cheap"),
                                        iterative=FALSE, precision=1e-3, Nmax=1e4, dN=Nquad,
                                        full=TRUE, progress=FALSE, verbose=TRUE,
                                        result.matrix=FALSE){

  averaging <- match.arg(averaging)

  wavelength <- material[["wavelength"]]
  k0 <- 2*pi/wavelength
  kn <- k0*medium
  
  Beta <- inverse_polarizability(cluster, material, 
                                 polarizability_fun=polarizability_ellipsoid, 
                                 medium=medium, kuwata=TRUE)
  
  
  quadrature <- integration_points(averaging, Nquad)
  
  results <- cd$average_spectrum(kn, Beta, cluster$r, cluster$angles, 
                                 as.matrix(quadrature$angles), 
                                 quadrature$weights,
                                 full, progress)
  
  ## iterative improvement: add new points until convergence or Nmax reached
  if(iterative && averaging == "QMC"){
    converged <- FALSE
    Ntot <- Nquad
    while(Ntot < Nmax && !converged){
      oldN <- Ntot
      old <- results[,1]
      Ntot <- Ntot + dN
      quadrature <- integration_points(averaging, dN, FALSE)
      ## xsec at new points
      newres <-  cd$average_spectrum(kn, Beta, cluster$r, cluster$angles, 
                                     as.matrix(quadrature$angles), 
                                     quadrature$weights,
                                     full, progress)
      
      ## average of the two results
      results <- (oldN * results + dN * newres) / (oldN + dN)
      
      test <- max(abs(old - results[,1]) / results[,1]) 
      ## max relative difference in extinction cross section
      if(verbose)
         message("N:", Ntot, "; relative error: " , test)
      converged <- test < precision
    }
  }
  
  
  if(result.matrix){
    ## extinction, absorption, scattering, CD ext, CD abs, CD sca
    return(cbind(wavelength = wavelength, results))
    
  } else {
    
    d <- data.frame(wavelength, results) # L - R
    names(d) <- c("wavelength", 'extinction', 'absorption', 'scattering',
                  "CDext", "CDabs", "CDsca")
    L2eV <- 6.62606896e-34 * 299792458/1.602176487e-19
    m <- melt(transform(d, energy = L2eV / wavelength * 1e9), id=c("wavelength", "energy"))
    
    m$type <- m$variable
    
    levels(m$type) <- list(CD="CDext",CD="CDabs",CD="CDsca",
                           `cross section`="extinction",
                           `cross section`="absorption",
                           `cross section`="scattering")
    
    levels(m$variable) <- list(extinction="extinction",
                               absorption="absorption",
                               scattering="scattering",
                               extinction="CDext",
                               absorption="CDabs",
                               scattering="CDsca")
    
    return(m)
    
  }

}

