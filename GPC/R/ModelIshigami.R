Model <- function(x){
  y=sin(2*pi*x[,1]-pi)+7*sin(2*pi*x[,2]-pi)^2+0.1*(2*pi*x[,3]-pi)^4*sin(2*pi*x[,1]-pi)
  return(y)
}