`foutliers` <-
function(data, method = c("robMah", "lrt", "depth.trim", "depth.pond", "HUoutliers"), 
             dfunc = depth.mode, nb = 200, suav = 0.05, trim = 0.1, order = 2, lambda = 3.29, ...)
{
  method = match.arg(method)
  if (method == "lrt"){
      output = outliers.lrt(data = data, nb = nb, suav = suav, trim = trim, ...)
  }
  if (method == "depth.trim"){
      output = outliers.depth.trim(data = data, dfunc = dfunc, nb = nb, suav = suav, trim = trim, ...)
  }
  if (method == "depth.pond"){
      output = outliers.depth.pond(data = data, dfunc = dfunc, nb = nb, suav = suav, ...)
  }
  if (method == "HUoutliers"){
      k = rapca(data$y, order = order)
      result = matrix(NA, ncol(data$y), 1)
      for(i in 1:ncol(data$y)){
          result[i,] = sum((data$y[,i] - k$basis %*% as.matrix(k$coef[i,]))^2)
      }
      s = median(result)
      crit = s+sqrt(s) * lambda
      out = as.numeric(colnames(data$y))[which(ifelse(result <= crit,1,0) == 0)]
      output = list(outliers = out)
  }
  if (method == "robMah"){
      sco = PCAproj(t(data$y))$scores
      rownames(sco) = as.numeric(colnames(data$y))
      s = cbind(sco, rep(1, ncol(data$y)))
      output = robout(s, 1, "mcd")
  }
  return(output)
}



