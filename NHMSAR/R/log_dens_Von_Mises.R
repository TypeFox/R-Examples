log_dens_Von_Mises <-
function (x, m, k) 
{   
	f = k*cos(x-m)-log(besselI(k,0))-log(2*pi)
}
