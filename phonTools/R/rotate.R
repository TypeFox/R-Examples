# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


rotate = function (xy, angle, degrees = FALSE, origin = TRUE){
  complx = FALSE
  if (degrees) angle = angle * pi/180
  if (is.complex(xy)) {
    xy = cbind(Re(xy), Im(xy))
    complx = TRUE
  }
  if (ncol(xy) != 2) stop("Input must have two columns (2-dimensional).")

  if (!origin){
    mus = colMeans (xy)
    xy = xy - matrix (mus,nrow(xy),2,byrow = T)
  }
  
  rotmat = matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
  output = as.matrix(xy) %*% rotmat
  
  if (!origin) output = output + matrix (mus,nrow(output),2,byrow = T)
  
  if (complx) output = complex(real = output[, 1], imaginary = output[,2])
  output
}

