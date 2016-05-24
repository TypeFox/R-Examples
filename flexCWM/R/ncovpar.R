ncovpar <- function (modelname = NULL, p = NULL, G = NULL) 
{
  if (is.null(p)) 
    stop("p is null")
  if (is.null(G)) 
    stop("G is null")
  if (is.null(modelname)) 
    stop("modelname is null")
  if (modelname == "EII") 
    npar = 1
  else if (modelname == "VII") 
    npar = G
  else if (modelname == "EEI") 
    npar = p
  else if (modelname == "VEI") 
    npar = p + G - 1
  else if (modelname == "EVI") 
    npar = p * G - G + 1
  else if (modelname == "VVI") 
    npar = p * G
  else if (modelname == "EEE") 
    npar = p * (p + 1)/2
  else if (modelname == "EEV") 
    npar = G * p * (p + 1)/2 - (G - 1) * p
  else if (modelname == "VEV") 
    npar = G * p * (p + 1)/2 - (G - 1) * (p - 1)
  else if (modelname == "VVV") 
    npar = G * p * (p + 1)/2
  else if (modelname == "EVE") 
    npar = p * (p + 1)/2 + (G - 1) * (p - 1)
  else if (modelname == "VVE") 
    npar = p * (p + 1)/2 + (G - 1) * p
  else if (modelname == "VEE") 
    npar = p * (p + 1)/2 + (G - 1)
  else if (modelname == "EVV") 
    npar = G * p * (p + 1)/2 - (G - 1)
  else stop("modelname is not correctly defined")
  return(npar)
}