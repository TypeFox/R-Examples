PROSPECT <- function(
                      N=1.5,
                      Cab=40,
                      Car=8,
                      Cbrown=0.0,
                      Cw=0.01,
                      Cm=0.009,
                      transmittance = FALSE,
                      parameterList = NULL
                   )
{
  if (!is.null(parameterList))
  {
    iterate_prospect <- function(x, transmittance = FALSE)
    {
      spec <- PROSPECT(N=x[1],
                       Cab=x[2],
                       Car=x[3],
                       Cbrown=x[4],
                       Cw=x[5],
                       Cm=x[6],
                       transmittance = transmittance)
      return(unlist(spectra(spec)[1,]))
    }
    
    parameter <- c("N", "Cab", "Car", "Cbrown", "Cw", "Cm")
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
    spec <- t(apply(parameterList, 1, FUN = iterate_prospect, 
                    transmittance = transmittance))
    return(speclib(spectra=spec, wavelength=c(1:2101)+399, 
                   attributes = as.data.frame(parameterList)))
  }

  nw=2101
  RT <- array(0, dim=c(nw,2))
  
  storage.mode(nw)     <- "integer" 
  storage.mode(N)      <- "double"
  storage.mode(Cab)    <- "double"
  storage.mode(Car)    <- "double"
  storage.mode(Cbrown) <- "double"
  storage.mode(Cw)     <- "double"
  storage.mode(Cm)     <- "double"
  storage.mode(RT)     <- "double"
  
  
  extern <- .Fortran("prospect2r",
                     N=N,
                     Cab=Cab,
                     Car=Car,
                     Cbrown=Cbrown,
                     Cw=Cw,
                     Cm=Cm,
                     RT2R=RT,
                     PACKAGE="hsdar"
             )
  
  spec <- speclib(wavelength = c(1:nw) + 399, 
                  spectra = matrix(data = extern$RT2R[, 1 + transmittance*1],
                                   nrow = 1)
                 )
  return(spec)
}
