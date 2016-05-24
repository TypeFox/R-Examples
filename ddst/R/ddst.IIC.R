`ddst.IIC` <-
function(coord, n, c=2.4) {
 IIC = as.numeric(max(coord[1],diff(coord)) < c*log(n))
 which.max(coord - (IIC*log(n) + (1-IIC)*2)*(1:length(coord))) 
}

