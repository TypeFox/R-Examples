#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


charStable <- function(theta, tau, pm = 0)
	{
	# pm is the type parametrization as described by Nolan(2009)
	# It takes the value 0 or 1 

	# const can fixe parameters. It is NULL for no constraint or
	# a matrix in which case the constraint is theta[const[,1]]=const[,2]

	a <- theta[1]
	b <- theta[2]
	g <- theta[3]
	d <- theta[4]
	if(pm == 0)
		{
		if(a == 1)
			{
			if(g == 0)
				{
				the_car <- exp(complex(imaginary = d*tau)) 
				}
			else
				{
				re_p <- -g * abs(tau)
				im_p <- d * tau
				im_p[tau!=0] <- im_p[tau != 0] + re_p[tau != 0]*2/pi*b*sign(tau[tau != 0])*log(g*abs(tau[tau != 0]))
				the_car <- exp(complex(real = re_p, imaginary = im_p))
				}
			}
		else
			{
			if(g == 0)
				{
				the_car <- exp(complex(imaginary = d*tau)) 
				}
			else
				{
				phi <- tan(pi*a/2)
				re_p <- -g^a*abs(tau)^a
				im_p <- d*tau*1i
				im_p[tau != 0] <- im_p[tau != 0] + re_p*( b*phi*sign(tau[tau != 0])*(abs(g*tau[tau != 0])^(1-a) - 1) )
				the_car <- exp(complex(real = re_p, imaginary = im_p))
				}
			}
		}

	if(pm == 1)
		{
		if(a == 1)
			{
			re_p <- -g*abs(tau)
			im_p <- d*tau
			im_p[tau!=0] <- im_p[tau != 0]+re_p*(b*2/pi*sign(tau[tau != 0])*log(abs(tau[tau!=0])))			
			the_car <- exp(complex(real = re_p, imaginary = im_p))
			}
		else
			{
			phi <- tan(pi*a/2)
			re_p <- -g^a*abs(tau)^a
			im_p <- re_p*(-phi*b*sign(tau)) + d*tau
			the_car <- exp(complex(real = re_p, imaginary = im_p))
			}
		}
	return(the_car)
	}



