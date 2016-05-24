# file eillks1998.R
#
# This file contains a script example with a precipitation stochastic generation 
#
#
# author: Emanuele Cordano on 10-02-2011

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################


library(RMRAINGEN)

p0_v <- c(0.69,0.73)
rho_range <- c(0,1)

x <- seq(from=min(rho_range),to=max(rho_range),by=0.001)

xrange <- range(x)
yrange <- xrange
xlab <- "Omega"
ylab <- "Xi"
main <- "Correlation dependence like Fig. 1 , Wilks, 1998, JoH"
cnt=0
for (p0_v1 in p0_v) {
	for (p0_v2 in p0_v) {
		cnt=cnt+1
		
		y <- omega(x=x,p0_v1=p0_v1,p0_v2=p0_v2,correlation=TRUE)
		if (cnt==1) {
			
			plot(x,y,xlim=xrange,ylim=yrange,type="l",lty=cnt,xlab=xlab,ylab=ylab,main=main)
		} else {
			lines(x,y,lty=cnt)
		}
	}
	
}
