#This function generates the specified number of points in the d-dim unit ball centred at the origin and radius 1
#For lower complexity, we now have a separate pathway for 2-dimensions
#Points are returned in polar coordinates
gen_points <- function(n,d)
{
	s <- sobol(n,dim=d,scrambling=3)
	r <- s[,1] ^ (1/d)	#The radius parameter
	if(d==2)
	{
		ang <- 2*pi*s[,2]
		return(cbind(r,ang))
	}
	else
	{
		ANGLES <- matrix(0,nrow=n,ncol=d-1)
		if(d>3)
		{
			construct_integrand <- function(k) g <- function(x) sin(x)^k
			k <- d - (1:(d-3)) - 1
			multiplier <- sqrt(pi) * gamma(k/ 2 + 0.5) / gamma(k/ 2 + 1)
			for(i in 1:(d-3))
			{
				integrand <- construct_integrand(k[i])
				rhs <- multiplier[i] * s[,i+1]
				for(j in 1:n)
				{
					f <- function(x) integrate(integrand,0,x)$value - rhs[j]
					ANGLES[j,i] <- uniroot(f , lower=0, upper = pi)$root
				}
			}
		}
		ANGLES[,d-2] <- acos(1 - 2*s[,d-1])
		ANGLES[,d-1] <- s[,d] * 2 * pi
	}
		return(cbind(r,ANGLES))
}
