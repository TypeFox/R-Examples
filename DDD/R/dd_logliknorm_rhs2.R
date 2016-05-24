dd_logliknorm_rhs2 = function(t,x,m)
{
   nx = sqrt(length(x))
   dim(x) = c(nx,nx)
   xx = matrix(0,nx+2,nx+2)
   xx[2:(nx+1),2:(nx+1)] = x
   dx = m[[1]] * xx[1:nx,2:(nx+1)] + m[[2]] * xx[3:(nx+2),2:(nx+1)] + m[[4]] * xx[2:(nx+1),1:nx] + m[[5]] * xx[2:(nx+1),3:(nx+2)] - (m[[3]] + m[[6]]) * xx[2:(nx+1),2:(nx+1)]
   dim(dx) = c(nx^2,1)
   return(list(dx))
}