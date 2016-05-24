 dd_logliknorm_rhs3 = function(t,x,m)
{
   nx = round((length(x))^(1/3))
   dim(x) = c(nx,nx,nx)     
   xx = array(0,dim = c(nx+2,nx+2,nx+2))
   xx[2:(nx+1),2:(nx+1),2:(nx+1)] = x
   dx = m[[1]] * xx[1:nx,2:(nx+1),2:(nx+1)] + m[[2]] * xx[3:(nx+2),2:(nx+1),2:(nx+1)] + m[[4]] * xx[2:(nx+1),1:nx,2:(nx+1)] + m[[5]] * xx[2:(nx+1),3:(nx+2),2:(nx+1)] + m[[7]] * xx[2:(nx+1),2:(nx+1),1:nx] + m[[8]] * xx[2:(nx+1),2:(nx+1),3:(nx+2)] - (m[[3]] + m[[6]] + m[[9]]) * xx[2:(nx+1),2:(nx+1),2:(nx+1)]
   dim(dx) = c(nx^3,1)
   return(list(dx))
}