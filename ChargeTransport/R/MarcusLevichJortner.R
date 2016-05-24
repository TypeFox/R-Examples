MarcusLevichJortner <- function(J, lambdaI, lambdaS, hBarW, dE0 = 0, dEField = 0, temp = 300, nVib = 50){
  if(missing(J))
    stop("Please specify the electronic coupling (J)")
  if(missing(lambdaI))
    stop("Please specify the internal reorganization energy (lambdaI)")
  if(missing(lambdaS))
    stop("Please specify the external reorganization energy (lambdaS)")
  if(missing(hBarW))
    stop("Please specify the energy of the effective vibrational mode (hBarW)")
    
  args <- list(
    J = J,
    lambdaI = lambdaI,
    lambdaS = lambdaS,
    hBarW = hBarW,
    dE0 = dE0,
    dEField = dEField,
    temp = temp,
    nVib = nVib)
  
  are.list <- sapply(args, is.list)
  if(any(are.list))
    stop("Arguments must be vectors, matrices or data.frames")
  
  are.data.frame <- sapply(args, is.data.frame)
  args[are.data.frame] <- lapply(args[are.data.frame], as.matrix)
  
  are.numeric <- sapply(args, is.numeric)
  if(any(!are.numeric))
    stop("Arguments must be numeric")
  
  are.scalar <- sapply(args, function(x) return(length(x) == 1))
  if(!are.scalar["temp"])
    stop("'temp' must be a scalar")
  if(!are.scalar["nVib"])
    stop("'nVib' must be a scalar")
  if(!are.scalar["hBarW"])
    stop("'hBarW' must be a scalar")
  
  lengths <- sapply(args[!are.scalar], length)
  dims    <- lapply(args[!are.scalar], dim   )
  
  if(length(lengths) != 0)
    if(!all(sapply(dims, identical, dims[[1]])) | !all(sapply(lengths, identical, lengths[[1]])))
      stop("Arguments must be scalars or have the same dimensions and lengths")
  
  hbar <- Joule2electronVolt(universalConstants["hbar","Value"]) # eV.s
  kb   <- Joule2electronVolt(universalConstants["kb"  ,"Value"]) # eV.K-1
  
  k <- 2*pi/hbar
  k <- k*(args$J)^2
  k <- k*sqrt(1/(4*pi*args$lambdaS*kb*args$temp))
  
  huangRhys = args$lambdaI/args$hBarW
  
  S <- 0
  for(n in 0:args$nVib)
  {
    toAdd <- exp(-huangRhys)
    toAdd <- toAdd*(huangRhys^n)/(factorial(n))
    toAdd <- toAdd*exp(-(args$dE0+args$dEField+args$lambdaS+args$hBarW*n)^2/
                         (4*args$lambdaS*kb*args$temp))
    S <- S + toAdd
  }
  k <- k*S
  
  attr(k, "rateExp") <- "Marcus-Levich-Jortner"
  
  return(k)
}