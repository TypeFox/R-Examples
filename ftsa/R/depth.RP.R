`depth.RP` <- function(data, trim = 0.25, nproj = 50, xeps = 0.0000001, x = NULL, graph = FALSE){
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
      mtrim[j,] = apply(functions[lista,],2,mean)
  }
  colnames(mtrim) = data$x
  mtrim = as.data.frame(mtrim, row.names = "")   
  if(graph){cgray = 1 - (prof - min(prof)) / (max(prof) - min(prof))
  if(p == 2){
      plot(range(functions[,1]), range(functions[,2]), type = "n")
      text(functions[,1], functions[,2], round(prof,3), cex = 0.75)
      points(med[1], med[2], col = 3, pch = 20, cex = 2)
      points(mtrim[,1], mtrim[,2], pch = 19, col = gray(2 * trim), cex = 2)
  } 
  else{
      plot(range(x), range(functions), type = "n", xlab = "t", ylab = "X(t)", main = "RP Depth")
      for(i in 1:n){
          lines(x, functions[i,], col = gray(cgray[i]))
      }
      lines(x, mtrim, lwd=2, col = "yellow")
      lines(x,med, col="red", lwd = 2) 
  }
  legend("topleft", legend = c(paste("Trim", trim*100, "\u0025", sep = ""), "Median"), lwd = 2, col =
         c("yellow", "red"))
  }  
  return(list("median" = med, "lmed" = k, "mtrim" = mtrim, "ltrim" = lista, "prof" = prof, "proj" = z))
}

