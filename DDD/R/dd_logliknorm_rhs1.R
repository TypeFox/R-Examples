dd_logliknorm_rhs1 = function(t,x,m)
{
   nx = length(x)
   m1 = m[1:nx]
   m2 = m[(nx+1):(2*nx)]
   m3 = m[(2*nx+1):(3*nx)]
   xx = c(0,x,0)
   dx = m1 * xx[1:nx] + m2 * xx[3:(nx+2)] - m3 * xx[2:(nx+1)]
   return(list(dx))
}