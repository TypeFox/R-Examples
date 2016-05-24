Marcus <- function(J, lambda, dE0 = 0, dEField = 0, temp = 300){
  if(missing(J))
    stop("Please specify the electronic coupling (J)")
  if(missing(lambda))
    stop("Please specify the reorganization energy (lambda)")
  
  args <- list(
    J = J,
    lambda = lambda,
    dE0 = dE0,
    dEField = dEField,
    temp = temp)
  
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

  lengths <- sapply(args[!are.scalar], length)
  dims    <- lapply(args[!are.scalar], dim   )

  if(length(lengths) != 0)
    if(!all(sapply(dims, identical, dims[[1]])) | !all(sapply(lengths, identical, lengths[[1]])))
      stop("Arguments must be scalars or have the same dimensions and lengths")
  
  hbar <- Joule2electronVolt(universalConstants["hbar","Value"]) # eV.s
  kb   <- Joule2electronVolt(universalConstants["kb"  ,"Value"]) # eV.K-1
  
  k <- 2*pi/hbar
  k <- k*(args$J)^2
  k <- k*sqrt(1/(4*pi*args$lambda*kb*args$temp))
  k <- k*exp(-(args$lambda+args$dE0+args$dEField)^2/(4*args$lambda*kb*args$temp))
  
  attr(k, "rateExp") <- "Marcus"
  
  return(k)
}