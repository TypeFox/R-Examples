##' permittivity silver
##'
##' analytical dielectric function of Silver (Drude model)
##' @title epsAg
##' @param wavelength wavelength in nm
##' @param epsilon.inf background dielectric constant
##' @param lambda.p plasma wavelength
##' @param mu.p damping constant
##' @return data.frame
##' @author baptiste Auguie
##' @export
##' @family user_level permittivity
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
##' @examples
##' require(dielectric) ; data(AgPalik)
##' wvl <- seq(300, 900)
##' silver <- epsAg(wvl)
##' 
##' matplot(silver$wavelength, cbind(Re(silver$epsilon), Im(silver$epsilon)), 
##' t="l", lty=1, xlab = "wavelength / nm", ylab = "Dielectric function")
##' matpoints(AgPalik$wavelength, cbind(Re(AgPalik$epsilon), Im(AgPalik$epsilon)), pch=1)

epsAg <- function(wavelength, epsilon.inf = 4,
                  lambda.p = 282, mu.p = 17000){
  
  data.frame(wavelength=wavelength, epsilon=
               epsilon.inf*(1 - 1 / (lambda.p^2*(1/wavelength^2 + 1i / (mu.p*wavelength)))))
}

##' permittivity gold
##'
##' analytical dielectric function of Au (Drude model + interband transitions)
##' @title epsAu
##' @param wavelength wavelength in nm
##' @param epsilon.infty background dielectric constant
##' @param lambda.p plasma wavelength
##' @param mu.p damping constant
##' @param A1 A1
##' @param phi1 phi1
##' @param lambda1 lambda1
##' @param mu1 mu1
##' @param A2 A2
##' @param phi2 phi2
##' @param lambda2 lambda2
##' @param mu2 mu2
##' @return data.frame
##' @export
##' @family user_level permittivity
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
##' @author baptiste Auguie
##' @examples
##' require(dielectric) ; data(AuJC)
##' wvl <- seq(300, 900)
##' gold <- epsAu(wvl)
##' 
##' matplot(gold$wavelength, cbind(Re(gold$epsilon), Im(gold$epsilon)), 
##' t="l", lty=1, xlab = "wavelength / nm", ylab = "Dielectric function")
##' matpoints(AuJC$wavelength, cbind(Re(AuJC$epsilon), Im(AuJC$epsilon)), pch=1)

epsAu <- function(wavelength, epsilon.infty = 1.54,
                  lambda.p = 177.5, mu.p = 14500,
                  A1 = 1.27, phi1 = -pi/4, lambda1 = 470, mu1 = 1900,
                  A2 = 1.1, phi2 = -pi/4, lambda2 = 325, mu2 = 1060){
  eps.drude <- 
    epsilon.infty*(1 - 1 / (lambda.p^2*(1/wavelength^2 + 1i / (mu.p*wavelength))))
  
  data.frame(wavelength=wavelength, epsilon=
               eps.drude + A1 / lambda1 * (exp(1i*phi1)  / (1/lambda1 - 1/wavelength - 1i/mu1) +
                                             exp(-1i*phi1) / (1/lambda1 + 1/wavelength + 1i/mu1)) +
               A2 / lambda2 * (exp(1i*phi2)  / (1/lambda2 - 1/wavelength - 1i/mu2) +
                                 exp(-1i*phi2) / (1/lambda2 + 1/wavelength + 1i/mu2))
  )
}
