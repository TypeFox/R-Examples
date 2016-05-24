`depth.RP` <- function(data, trim = 0.25, nproj = 50, xeps = 0.0000001, x = NULL){
  functions = t(data$y)
  n = dim(functions)[1]
  p = dim(functions)[2]
  if(is.null(x)) 
     x = 1:p
  prof = rep(0.0,n)
  for(j in 1:nproj){
      z = rnorm(p)
      modulo = sum(z^2)
      z = z/sqrt(modulo)
      valor = functions %*% z
      Fn = ecdf(valor)
      prof = prof + (Fn(valor) * (1 - Fn(valor - xeps)))
  }
  prof = prof / nproj
  k = which.max(prof)
  med = functions[k,]
  nl = length(trim)
  mtrim = matrix(NA, nrow = nl, ncol = dim(functions)[2])
  for(j in 1:length(trim)) {
      lista = which(prof >= quantile(prof, probs = trim[j], na.rm = TRUE))
      mtrim[j,] = apply(functions[lista,], 2, mean)
  }
  return(list("median" = med, "lmed" = k, "mtrim" = mtrim, 
         "ltrim" = lista, "prof" = prof))
}

