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

data(trentino) 

year_min <- 1961
year_max <- 1990

period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max

###period <- period & (PRECIPITATION$month %in% c(9,10,11)) ## AUTUMN 

station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
prec_mes <- PRECIPITATION[period,station]

## removing nonworking stations (e.g. time series with NA)
accepted <- array(TRUE,length(names(prec_mes)))
names(accepted) <- names(prec_mes)
for (it in names(prec_mes)) {
	accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
}

prec_mes <- prec_mes[,accepted]
names <- c("T0083","T0090","T0129") ###"T0083",
prec_mes <- prec_mes[,names(prec_mes) %in% names]
names <- names(prec_mes)



rho_range <- c(0,1)

x <- seq(from=min(rho_range),to=max(rho_range),by=0.01)

xrange <- range(x)
yrange <- xrange
xlab <- "Sigma_z"
ylab <- "Xi"
main <- "Xi Correlation vs Gaussian z Correlation"#"Correlation dependence like Fig. 1 , Wilks, 1998, JoH"
cnt=1
legend <- NA
type <- NA

seasons <- list(DJF=c(12,1,2),MAM=c(3,4,5),JJA=c(6,7,8),SON=c(9,10,11)) 
for (seas in names(seasons)) {	
	
	season <- PRECIPITATION$month %in% seasons[[seas]] 
	
	p0_v <- diag(continuity_ratio(prec_mes[season,],lag=0)$nooccurence)
	
	names(p0_v) <-  names
	
	print(seas)
	print(p0_v)
	for (l in names) { 
		
		
		p0_v1 <- p0_v[l]
		
		for (c in names) { 
			
			p0_v2 <- p0_v[c] 
			
			
			y <- omega(x=x,p0_v1=p0_v1,p0_v2=p0_v2,correlation=TRUE)
			if ((cnt==1) & (p0_v1<p0_v2)) {
				
				plot(x,y,xlim=xrange,ylim=yrange,type="l",lty=cnt,xlab=xlab,ylab=ylab,main=main)
				legend[cnt] <- paste(l,c,sep="-")
				legend[cnt] <- paste(legend[cnt]," (",seas,")",sep="")
				type[cnt] <- cnt
				cnt=cnt+1
			} else if (p0_v1<p0_v2) {
				lines(x,y,lty=cnt)
				legend[cnt] <- paste(l,c,sep="-")
				legend[cnt] <- paste(legend[cnt]," (",seas,")",sep="")
				type[cnt] <- cnt
				cnt=cnt+1
			}
		}
		
	}
}
legend("topleft",lty=type,legend=legend)