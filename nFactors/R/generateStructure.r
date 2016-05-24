generateStructure <-
function(var, mjc, pmjc, loadings, unique) {

 if (var/mjc != ceiling(var/mjc)) stop("Bad pmjc value")
 
 fload <- matrix(unique, ncol=mjc, nrow=var)
 
 for (i in 1:mjc) {
  if (i == 1) fload[i:(pmjc),i] <- loadings
  min <- ((i-1)*pmjc+1)
  max <- ((i-1)*pmjc+pmjc)
  if (min > dim(fload)[1]) min <- dim(fload)[1]
  if (max > dim(fload)[1]) max <- dim(fload)[1]
  if (i > 1)  fload[min:max,i] <- loadings
  if (min > dim(fload)[1]) fload[dim(fload)[1],i] <- unique
  }
  
 return(data.frame(fload))
 }