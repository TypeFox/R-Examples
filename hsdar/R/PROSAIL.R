PROSAIL <- function(
                    N=1.5,
                    Cab=40,
                    Car=8,
                    Cbrown=0.0,
                    Cw=0.01,
                    Cm=0.009,
                    psoil=0,
                    LAI=1,
                    TypeLidf= 1,
                    lidfa= -0.35,
                    lidfb= -0.15,
                    hspot=0.01,
                    tts = 30,
                    tto = 10, 
                    psi= 0,
                    parameterList = NULL,
                    rsoil = NULL
                   )
{
  if (!is.null(parameterList))
  {
    iterate_prosail <- function(x, rsoil)
    {
      spec <- PROSAIL(N        = x[1],
                      Cab      = x[2],
                      Car      = x[3],
                      Cbrown   = x[4],
                      Cw       = x[5],
                      Cm       = x[6],
                      psoil    = x[7],
                      LAI      = x[8],
                      TypeLidf = x[9],
                      lidfa    = x[10],
                      lidfb    = x[11],
                      hspot    = x[12],
                      tts      = x[13],
                      tto      = x[14], 
                      psi      = x[15],
                      rsoil    = rsoil
                     )
      return(unlist(spectra(spec)[1,]))
    }
    
    parameter <- c("N", "Cab", "Car", "Cbrown", "Cw", "Cm",
                   "psoil", "LAI", "TypeLidf", "lidfa", 
                   "lidfb", "hspot", "tts", "tto", "psi")
    parameterList <- as.data.frame(parameterList)
    nam_para <- names(parameterList)
    mat <- match(names(parameterList), parameter, nomatch=0)
    if (any(mat==0))
      stop("Check names and format of parameterList")
    mat <- match(c(1:length(parameter)), mat, nomatch=0)
    for (i in 1:length(mat))
      if (mat[i]==0)
        parameterList <- cbind(parameterList, get(eval(parameter[i])))
    names(parameterList) <- c(nam_para, parameter[mat==0])
    
    parameterList <- as.matrix(parameterList[,match(parameter, names(parameterList))])
    spec <- t(apply(parameterList, 1, FUN = iterate_prosail, rsoil))
    return(speclib(spectra=spec, wavelength=c(1:2101)+399, 
                   attributes = as.data.frame(parameterList)))
  }
  nw <- 2101
  RT <- array(0, dim=c(nw,1))
  
  if (!is.null(rsoil))
  {
    if (length(wavelength(rsoil)) != nw)
      stop("Wavelength of rsoil must be 400 - 2500 nm with 1 nm resolution")
    if (!all(wavelength(rsoil) == c(400:2500)))
      stop("Wavelength of rsoil must be 400 - 2500 nm with 1 nm resolution")
    if (nspectra(rsoil) > 1)
      warning("More than one spectrum in Speclib rsoil. Only the first one will be used for simulation")
    rsoil <- spectra(rsoil)[1,]
  } else {
    rsoil <- rep.int(-9999.9, nw)
  }
  
  if (TypeLidf != 1) TypeLidf <- 2
    
  if (TypeLidf == 2) lidfb <- 0
  
  storage.mode(nw)       <- "integer" 
  storage.mode(N)        <- "double"
  storage.mode(Cab)      <- "double"
  storage.mode(Car)      <- "double"
  storage.mode(Cbrown)   <- "double"
  storage.mode(Cw)       <- "double"
  storage.mode(Cm)       <- "double"
  storage.mode(psoil)    <- "double"
  storage.mode(LAI)      <- "double"
  storage.mode(TypeLidf) <- "integer" 
  storage.mode(lidfa)    <- "double"
  storage.mode(lidfb)    <- "double"
  storage.mode(hspot)    <- "double"
  storage.mode(tts)      <- "double"
  storage.mode(tto)      <- "double"
  storage.mode(psi)      <- "double"
  storage.mode(RT)       <- "double"  
  storage.mode(rsoil)    <- "double"  
  
  extern <- .Fortran("prosail2r",
                     Cab=Cab,
                     Car=Car,
                     Cbrown=Cbrown,
                     Cw=Cw,
                     Cm=Cm,
                     N=N,
                     psoil=psoil,
                     LAI=LAI,
                     TypeLidf=TypeLidf,
                     lidfa=lidfa,
                     lidfb=lidfb,
                     hspot=hspot,
                     tts=tts,
                     tto=tto,
                     psi=psi,
                     reflectance=RT,
                     rsoil=rsoil,
                     PACKAGE="hsdar"
             )
  spec <- speclib(wavelength=c(1:nw)+399,
                  spectra=extern$reflectance
                 )
  return(spec)
}
