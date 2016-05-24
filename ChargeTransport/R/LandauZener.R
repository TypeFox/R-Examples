LandauZener <- function(J, lambda, nuN, dE0 = 0, dEField = 0, temp = 300, gammaN = 1){
  if(missing(J))
    stop("Please specify the electronic coupling (J)")
  if(missing(lambda))
    stop("Please specify the reorganization energy (lambda)")
  if(missing(nuN))
    stop("Please specify the frequency of the effective vibrational mode (nuN)")
  
  args <- list(
    J = J,
    lambda = lambda,
    nuN = nuN,
    dE0 = dE0,
    dEField = dEField,
    temp = temp,
    gammaN = gammaN)
  
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
  if(!are.scalar["nuN"])
    stop("'nuN' must be a scalar")
  if(!are.scalar["gammaN"])
    stop("'gammaN' must be a scalar")
  
  lengths <- sapply(args[!are.scalar], length)
  dims    <- lapply(args[!are.scalar], dim   )
  
  if(length(lengths) != 0)
    if(!all(sapply(dims, identical, dims[[1]])) | !all(sapply(lengths, identical, lengths[[1]])))
      stop("Arguments must be scalars or have the same dimensions and lengths")
  
  kb   <- Joule2electronVolt(universalConstants["kb"  ,"Value"]) # eV.K-1
  
  kel   <- electronicTransmissionCoefficient(args$J, args$lambda, args$nuN, args$temp)
  dE    <- activationFreeEnergy(args$lambda, args$dE0 + args$dEField)
  Delta <- adiabaticCorrection(args$J, args$lambda, args$dE0 + args$dEField)
 
  k <- kel*args$nuN
  k <- k*args$gammaN
  k <- k*exp(-(dE-Delta)/(kb*args$temp))
 
  attr(k, "rateExp") <- "Landau-Zener"
  
  return(k)
}