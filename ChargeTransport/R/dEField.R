dEField <- function( dx, dy, dz, carrier = "e", Fx = 0, Fy = 0, Fz = 1E5)
{
  if(missing(dx)) stop("Please specify the x-component of the intersite distance")
  if(missing(dy)) stop("Please specify the y-component of the intersite distance")
  if(missing(dz)) stop("Please specify the z-component of the intersite distance")
  
  args <- list(
    dx = dx,
    dy = dy,
    dz = dz,
    Fx = Fx,
    Fy = Fy,
    Fz = Fz)
  
  are.list <- sapply(args, is.list)
  if(any(are.list))
    stop("Arguments must be vectors, matrices or data.frames")
  
  are.data.frame <- sapply(args, is.data.frame)
  args[are.data.frame] <- lapply(args[are.data.frame], as.matrix)
  
  are.numeric <- sapply(args, is.numeric)
  if(any(!are.numeric))
    stop("Arguments must be numeric")
  
  are.scalar <- sapply(args, function(x) return(length(x) == 1))
  if(!are.scalar["Fx"])
    stop("'Fx' must be a scalar")
  if(!are.scalar["Fx"])
    stop("'Fy' must be a scalar")
  if(!are.scalar["Fx"])
    stop("'Fz' must be a scalar")
  
  lengths <- sapply(args[!are.scalar], length)
  dims    <- lapply(args[!are.scalar], dim   )
  
  if(length(lengths) != 0)
    if(!all(sapply(dims, identical, dims[[1]])) | !all(sapply(lengths, identical, lengths[[1]])))
      stop("Arguments must be scalars or have the same dimensions and lengths")
    
  if(carrier=="h")
    siteEnergyDiff <- -(dx*Fx+dy*Fy+dz*Fz)*1E-8
  else if(carrier=="e")
    siteEnergyDiff <-  (dx*Fx+dy*Fy+dz*Fz)*1E-8
  else stop("Unrecognized type of charge carrier. 'carrier' must be equal to 'e' or 'h'")
  
  return(siteEnergyDiff)
}
