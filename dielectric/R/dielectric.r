
dielectric <- setRefClass("dielectric",
                          fields = list( wavelength = "vector",
                            epsilon = "vector",
                            span="vector", 
                            comment = "list"),
                          methods = list(
                            set_span = function(min=NULL, max=NULL){
                              if(is.null(min)) min <- min(wavelength)
                              if(is.null(max)) max <- max(wavelength)
                              span <<- c(min, max)
                            }, 
                            ##' spline interpolation of permittivity
                            ##'
                            ##' spline interpolation of permittivity
                            ##' @title spline
                            ##' @param ... 
                            ##' @return list
                            ##' @author baptiste Auguie
                            spline = function(...) {
                              'returns a list of splinefun for the real and imaginary parts'
                              sp <- 
                                list(real=smooth.spline(wavelength, Re(epsilon), ...),
                                     imag=smooth.spline(wavelength, Im(epsilon), ...))
                              invisible(sp)
                            },
                            raw = function(range=span){
                              'return the raw data as real numbers'
                            subset(data.frame(wavelength=wavelength,
                                       epsilon=epsilon, real = Re(epsilon), imag=Im(epsilon)),
                                   wavelength < range[2] & wavelength > range[1])
                          },
                            permittivity = function(new.wavelength, ...){
                              'predict a single value'
                              predict(range=rep(new.wavelength, 2), n=1, ...)
                              },
                            predict = function(sp=NULL, range=span,
                              n=length(epsilon), new.wavelength=NULL,...)
                            {
                              'interpolation with splines'
                              if(is.null(range))
                                 range <- range(wavelength)
                              
                              if(is.null(new.wavelength)){ 
                                new.wavelength <- if(is.null(n)) { # use the original points
                                  wavelength[findInterval(wavelength, range)] } else
                                {
                                  seq(range[1], range[2], length=n)
                                }
                              }
                              
                              sp <- if(is.null(sp)) spline(...)
                              
                              smooth <-
                                data.frame(wavelength=new.wavelength,
                                           epsilon = complex(real=stats::predict(sp$real,
                                                               new.wavelength)$y,
                                             imag=stats::predict(sp$imag, new.wavelength)$y))
                              invisible(smooth)
                            }
                            ))

