##' Multilayer Fresnel coefficients
##'
##' computes the reflection coefficient of a multilayered interface
##' @title recursive_fresnel
##' @export
##' @param wavelength [vector] wavelength in nm
##' @param k0 [vector] wavevector in nm^-1
##' @param angle [vector] incident angles in radians
##' @param q [vector] normalised incident in-plane wavevector
##' @param epsilon list of N+2 dielectric functions, each of length 1 or length(wavelength)
##' @param thickness vector of N+2 layer thicknesses, first and last are dummy
##' @param polarisation [character] switch between p- and s- polarisation
##' @return fresnel coefficients and field profiles
##' @author baptiste Auguie
recursive_fresnel <- function(wavelength = 2*pi/k0, k0 = 2*pi/wavelength,
                       angle = NULL, q = sin(angle),
                       epsilon = list(incident=1.5^2, 1.33^2),
                       thickness = c(0, 0),
                       polarisation = c('p', 's')){

  if(is.null(angle))
    suppressWarnings(angle <- asin(q))
  ## checks
  # stopifnot(thickness[1]==0L, thickness[length(thickness)]==0L)
  polarisation <- match.arg(polarisation)

  ## define constants
  Nlambda <- length(k0)
  Nq <- length(q)
  Nlayer <- length(thickness)
  k02 <- k0^2
  kx <- outer(k0*sqrt(epsilon[[1]] + 0i), q) # kx = q*k0
  kx2 <- kx^2

  ## loop to calculate kiz
  kiz <- array(0 + 0i, dim=c(Nlambda, Nq, Nlayer))

  for (ii in seq(1, Nlayer)){
    kiz[ , , ii] <- sqrt(matrix(epsilon[[ii]]*k02, nrow=Nlambda, ncol=Nq) - kx2 + 0i)
  }

  ## calculate all single-interface fresnel coefficients and phase factors

  Nsurface <- Nlayer -1
  rsingle <- tsingle <- array(0 + 0i, dim=c(Nlambda, Nq, Nsurface))
  phase1 <- phase2 <- array(1 + 0i, dim=c(Nlambda, Nq, Nlayer))

  for (ii in seq(1, Nsurface)){

    if(polarisation == 'p'){
      # Note: this is for the H field
      a <- kiz[,,ii] / epsilon[[ii]]
      b <- kiz[,,ii+1] / epsilon[[ii+1]]

    } else { # s-polarisation
      # for the E field
      a <- kiz[,,ii]
      b <- kiz[,,ii+1]

    }

    rsingle[,,ii] <-  (a - b) / (a + b)
    tsingle[,,ii] <-  2 * a / (a + b)

    phase1[,,ii] <- exp(1i*thickness[ii]*kiz[,,ii])
    phase2[,,ii] <- exp(2i*thickness[ii]*kiz[,,ii])

  }

  phase1[,,Nlayer] <- 1 + 0i # 0 thickness for last medium
  phase2[,,Nlayer] <- 1 + 0i # 0 thickness for last medium

  ## now recursion, r.tmp is the combined reflection from last layer
  ## r_{N-1,N}, then r_{N-2,N}, ... finally r_{1N}

  ## starting from last rsingle
  reflection <- rsingle[,,Nsurface]
  transmission <- tsingle[,,Nsurface]

  for(jj in rev(seq_len(Nsurface - 1))){
    transmission <- ( tsingle[,,jj] * transmission*phase1[,,jj+1]) /
      (1 + rsingle[,,jj]*reflection*phase2[,,jj+1])
    reflection <- ( rsingle[,,jj] + reflection*phase2[,,jj+1]) /
      (1 + rsingle[,,jj]*reflection*phase2[,,jj+1])
  }

  ## T is nt*cos(Ot)*|Et|^2 / ni*cos(Oi)*|Ei|^2
  ## for s-pol, |Et|^2 / |Ei|^2 = |ts|^2, hence T = nt/ni * cos(Ot)/cos(Oi) * |ts|^2
  ## for p-pol, |Et|^2 / |Ei|^2 = (ni/nt)^2 * |tp|^2, hence T = ni/nt * cos(Ot)/cos(Oi) * |tp|^2

  # ratio of refractive indices
  index.ratio <- matrix(Re(sqrt(epsilon[[1]])/sqrt(epsilon[[Nlayer]])), nrow=Nlambda, ncol=Nq)
  # ratio of cosines
  qq <- matrix(q, nrow=Nlambda, ncol=Nq, byrow=TRUE)
  m <- Re(sqrt(1 - (index.ratio * qq)^2 + 0i)/sqrt(1 - qq^2 + 0i))

  if(polarisation == "p"){
    rho <- index.ratio
    reflection <- -reflection # sign convention was for H
  } else {
    rho <- 1 / index.ratio
  }
  R <- Mod(reflection)^2
  T <- rho * m * Mod(transmission)^2

  list(wavelength=wavelength, k0=k0,
       angle=angle, q=q,
       reflection=reflection,
       transmission=transmission,
       R = R, T = T, A = 1 - R - T)

}

##' Multilayer Fresnel coefficients
##'
##' computes the reflection coefficient of a multilayered interface
##' @title recursive_fresnelcpp
##' @export
##' @param wavelength [vector] wavelength in nm
##' @param k0 [vector] wavevector in nm^-1
##' @param angle [vector] incident angles in radians
##' @param q [vector] normalised incident in-plane wavevector
##' @param epsilon list of N+2 dielectric functions, each of length 1 or length(wavelength)
##' @param thickness vector of N+2 layer thicknesses, first and last are dummy
##' @param polarisation [character] switch between p- and s- polarisation
##' @return fresnel coefficients and field profiles
##' @author baptiste Auguie
recursive_fresnelcpp <- function(wavelength = 2*pi/k0, k0 = 2*pi/wavelength,
                       angle = NULL, q = sin(angle),
                       epsilon = list(incident=1.5^2, 1.33^2),
                       thickness = c(0, 0),
                       polarisation = c('p', 's')){

  if(is.null(angle))
    suppressWarnings(angle <- asin(q))

  polarisation <- match.arg(polarisation)
  kx <- outer(k0*sqrt(epsilon[[1]]), q) # kx = q*k0
  epsilon = do.call(cbind, epsilon)

  Nlambda <- length(k0)
  Nq <- length(q)
  Nlayer <- length(thickness)

  ## case pure scalars
  if(nrow(epsilon) == 1L)
    epsilon <- matrix(epsilon, nrow=length(k0), ncol=Nlayer, byrow=TRUE)

  polarisation = if(polarisation == "p") 0L else 1L

  ## checks
  stopifnot(thickness[1]==0L,
            thickness[Nlayer]==0L)


  stopifnot(Nlayer == ncol(epsilon),
            nrow(epsilon) == length(k0),
            nrow(kx) == length(k0),
            ncol(kx) == length(q))


  ## call the C++ function
  res <- cpp_recursive_fresnel(as.vector(k0),
                                  as.matrix(kx),
                                  as.matrix(epsilon),
                                  as.vector(thickness),
                                  as.integer(polarisation))

  transmission <- drop(res$transmission)
  reflection <- drop(res$reflection)

  ## T is nt*cos(Ot)*|Et|^2 / ni*cos(Oi)*|Ei|^2
  ## for s-pol, |Et|^2 / |Ei|^2 = |ts|^2, hence T = nt/ni * cos(Ot)/cos(Oi) * |ts|^2
  ## for p-pol, |Et|^2 / |Ei|^2 = (ni/nt)^2 * |tp|^2, hence T = ni/nt * cos(Ot)/cos(Oi) * |tp|^2

  # ratio of refractive indices
  index.ratio <- matrix(Re(sqrt(epsilon[, 1])/sqrt(epsilon[, Nlayer])), nrow=Nlambda, ncol=Nq)
  # ratio of cosines
  qq <- matrix(q, nrow=Nlambda, ncol=Nq, byrow=TRUE)
  m <- Re(sqrt(1 - (index.ratio * qq)^2 + 0i)/sqrt(1 - qq^2 + 0i))


  if(polarisation == 0L){
    rho <- index.ratio
    reflection <- -reflection # sign convention was for H
  } else {
    rho <- 1 / index.ratio
  }
  R <- Mod(reflection)^2
  T <- rho * m * Mod(transmission)^2

  list(wavelength = wavelength, k0 = k0,
       angle=angle, q=q,
       reflection=reflection,
       transmission=transmission,
       R=R, T=T, A = 1 - R - T)
}



single_layer <- function(wavelength = 2*pi/k0, k0 = 2*pi/wavelength,
                         angle = NULL, q = sin(angle),
                         epsilon = list(incident=3.0^2, 3.7^2, 3.0^2),
                         thickness = 400){

  if(is.null(angle))
    suppressWarnings(angle <- asin(q))

  kx <- outer(k0*sqrt(epsilon[[1]]), q) # kx = q*k0
  epsilon = do.call(cbind, epsilon)

  Nlayer <- 3

  ## case pure scalars
  if(nrow(epsilon) == 1L)
    epsilon <- matrix(epsilon, nrow=length(k0), ncol=Nlayer, byrow=TRUE)

  stopifnot(Nlayer == ncol(epsilon),
            nrow(epsilon) == length(k0),
            nrow(kx) == length(k0),
            ncol(kx) == length(q))

  ## call the C++ function
  res <- cpp_layer_fresnel(as.vector(k0),
                              as.matrix(kx),
                              as.matrix(epsilon),
                              as.vector(thickness))
  res[['rp']]
}
