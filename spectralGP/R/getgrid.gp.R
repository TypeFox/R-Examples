"getgrid.gp" <-
function(object,...){
  # return the implicit grid of the process as two columns, the first the grid values in the x-direction and the second in the y-direction
  m1=object$gridsize[1]
  m2=object$gridsize[2]
  if(object$d==2){
    return(list(x1=seq(0,1,length=m1/2+1),x2=seq(0,1,length=m2/2+1)))
  } else{
    return(seq(0,1,length=m1/2+1))
  }
}
