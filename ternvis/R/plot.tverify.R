plot.tverify <-
function(
                       x,            # of class tverify
		               thresh=0,       # screen out bins with < thresh data points
		               lsharp=TRUE,    # switch to show sharpness diagram in top right
                       L = diag(c(1,1,1))/sqrt(2),...
		       )		       		           
{
#
# define utility functions
#
dat <- x
 
tsharp <- function(dat,          # object of class tverify
		      main="",
		      sub="",
		      lfigs=FALSE   # logical to print numbers inside hexagons
		      )
{

q <- dat$q

temp <- tsetup()
xlim <- temp$xlim
ylim <- temp$ylim
xoff <- temp$xoff
yoff <- temp$yoff
rm(temp)


# a nice colour palette from RColorBrewer:
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")

hexc   <- dat$hexc
dens   <- dat$Nobs
parray <- dat$pbin
bigN   <- (dat$ncirc+1)*(dat$ncirc+2)/2
cols   <- rep(NA,bigN)
L      <- dat$L

maxdens=max(dens) 
mindens=max(min(dens),1e-6) # a hack

for (i in 1:bigN){
   if (dens[i] == 0) 
      {cols[i] <- "lightgrey"}
   else
      {rgbtemp <- colorRamp(YlOrBr)((log(dens[i])-log(mindens))/(log(maxdens)-log(mindens)))/255
       cols[i] <- rgb(rgbtemp[1],rgbtemp[2],rgbtemp[3])}
}

x <- cbind(xf(hexc[,1,]),
           xf(hexc[,2,]),
	   xf(hexc[,3,]),
	   xf(hexc[,4,]),
	   xf(hexc[,5,]),
	   xf(hexc[,6,]))

y <- cbind(yf(hexc[,1,]),
           yf(hexc[,2,]),
	   yf(hexc[,3,]),
	   yf(hexc[,4,]),
	   yf(hexc[,5,]),
	   yf(hexc[,6,]))	   

tplot(q,
        col="transparent",
        L=L,
        dimnames=NULL,
        main=main,
	sub=sub,
	newpage=FALSE,
	grid=FALSE
	)
      
pushViewport(viewport(width  = 0.8, 
                      height = 0.8, 
		      xscale = xlim,
                      yscale = ylim))

j <- order(dens) # we will want to fill in polygons in acsending order of density

for (i in 1:bigN)
{
   grid.polygon(x=xoff+x[j[i],],
                y=yoff+y[j[i],],
		default.units="native",
		gp=gpar(col=cols[j[i]],fill=cols[j[i]]))
}

if (lfigs)
{
grid.text(paste(dens),
          x=xoff+xf(parray),
	  y=yoff+yf(parray),
	  default.units="native",
	  gp=gpar(cex=11/dat$ncirc))
}
popViewport()

tplot( q,
       L=L,
       col="blue", 
       pch=4, 
       cex=0.5,
       main=paste(""),
       newpage=FALSE,
       bg="transparent",
       grid=FALSE,
       dimnames_position="none",
       border="black")  # plot climatology
} # end of function tsharp


tdecomp <- function(
           dat,         # an object of class tverify
		   ldec=FALSE   # logical to show decimal numerical values of quantities, otherwise give symbol
		   )
{
uncbar <- dat$uncbar
resbar <- dat$resbar
relbar <- dat$relbar

temp   <- tsetup()
maxhw  <- temp$maxhw
width  <- temp$width
ymin   <- temp$ymin
height <- temp$height
yoff   <- temp$yoff
rm(temp)


scorebar = uncbar - resbar + relbar
size <- (maxhw/11)/5
xoff <- 0.0
toff <- width/30
ytop <- ymin+height*1.03
x1   <- 0                         
y1   <- ytop
x4   <- 0
y4   <- y1 - sqrt(uncbar)

theta1 <- asin(sqrt(resbar)/sqrt(uncbar))
theta2 <- asin(sqrt(relbar)/sqrt(scorebar))

theta <- theta1-theta2

x2 <- 0  + sqrt(uncbar-resbar) * sin(theta1)
y2 <- y4 + sqrt(uncbar-resbar) * cos(theta1)

x3 <- 0  + sqrt(scorebar) * sin(theta)
y3 <- y4 + sqrt(scorebar) * cos(theta)


        grid.polygon(x=xoff+c(
	                      x2,
	                      x2-toff*cos(theta1),
	                      x2-toff*cos(theta1)-toff*sin(theta1),
			      x2-toff*sin(theta1)
			      ),
                     y=yoff+c(
		              y2,
		              y2+toff*sin(theta1),
			      y2+toff*sin(theta1)-toff*cos(theta1),
			      y2-toff*cos(theta1)
			      ),
		     default.units="native",
		     gp=gpar(col="black"))



th <- seq(0,theta1,length.out=100)
xarc1  <- x4 + sqrt(uncbar)*sin(th)
yarc1  <- y4 + sqrt(uncbar)*cos(th)
xarc2  <- x4 + sqrt(uncbar-resbar)*sin(th)
yarc2  <- y4 + sqrt(uncbar-resbar)*cos(th)

        grid.lines(x=xoff+xarc1,
                   y=yoff+yarc1,
		   default.units="native",
		   gp=gpar(lwd=1,lty=2,col="blue"))        # arc of maximum score

        grid.lines(x=xoff+xarc2,
                   y=yoff+yarc2,
		   default.units="native",
		   gp=gpar(lwd=1,lty=2,col="purple"))      # arc of minimum score

th <- seq(0,pi,length.out=100)
xarc3  <- 0.5*(x1+x4) + 0.5 * sqrt(uncbar)*sin(th)
yarc3  <- 0.5*(y1+y4) + 0.5 * sqrt(uncbar)*cos(th)

        grid.lines(x=xoff+xarc3,
                   y=yoff+yarc3,
		   default.units="native",
		   gp=gpar(col="grey"))                     # semicircle








        grid.lines(x=xoff+c(x1,x2),
                   y=yoff+c(y1,y2),
		   default.units="native",
		   gp=gpar(col="darkgreen",lwd=1))

		        
        grid.lines(x=xoff+c(x2,x3),
                   y=yoff+c(y2,y3),
		   default.units="native",
		   gp=gpar(col="red",lwd=1))

        grid.lines(x=xoff+c(x1,x4),
                   y=yoff+c(y1,y4),
		   default.units="native",
		   gp=gpar(col="blue",lwd=1))

        grid.lines(x=xoff+c(x2,x4),
                   y=yoff+c(y2,y4),
		   default.units="native",
		   gp=gpar(col="purple",lwd=1))

		   
        grid.lines(x=xoff+c(x3,x4),
                   y=yoff+c(y3,y4),
		   default.units="native",
		   gp=gpar(col="black",lwd=3))
		   
		
	grid.points(x=xoff+x3,
	            y=yoff+y3,
		    pch=21,   # circle
		    size=unit(size,"native"),
		    default.units="native",
		    gp=gpar(col="black",fill="transparent"))	# end of score line

	grid.points(x=xoff+x2,
	            y=yoff+y2,
		    pch=21,   # circle
		    size=unit(size,"native"),		    
		    default.units="native",
		    gp=gpar(col="red",fill="transparent"))      # right angle point



toff <- width / 30 # an offset for the text labels
 
        label <- expression(sqrt(S))
	if (ldec) label <- sprintf("%.3f",sqrt((x3-x4)^2+(y3-y4)^2))
        grid.text(label=label,
	      x=xoff+0.7*x3+0.3*x4+toff*cos(theta),
		  y=yoff+0.7*y3+0.3*y4-toff*sin(theta),
		  default.units="native",
		  rot=90-(180/pi)*theta,
		  just="centre")


        label <- expression(sqrt(R))
        if (ldec) label <- sprintf("%.3f",sqrt((x2-x3)^2+(y2-y3)^2))
        grid.text(label=label,
	          x=xoff+x1 + sqrt(uncbar) * sin(theta + theta2/2.) + toff*sin(theta1),
		  y=yoff+y4 + sqrt(uncbar) * cos(theta + theta2/2.) + toff*cos(theta1),
		  default.units="native",
		  just="centre",
		  rot=-(180/pi)*(theta1),		  
		  gp=gpar(col="red"))

        label <- expression(sqrt(Z))
        if (ldec) label <- sprintf("%.3f",sqrt((x1-x2)^2+(y1-y2)^2))
        grid.text(label=label,
	          x=xoff+x1 + sqrt(uncbar) * sin(theta1/2.) + 3*toff*sin(theta1),
		  y=yoff+y4 + sqrt(uncbar) * cos(theta1/2.) + 3*toff*cos(theta1),
		  default.units="native",
		  just="centre",
		  rot=-(180/pi)*theta1,		  
		  gp=gpar(col="darkgreen"))

        label <- expression(sqrt(U))
        if (ldec) label <- sprintf("%.3f",sqrt((x1-x4)^2+(y1-y4)^2))
        grid.text(label=label,
	          x=xoff+x1-toff,
		  y=yoff+0.5*(y1+y4),
		  default.units="native",
		  just="centre",
		  rot=90,
		  gp=gpar(col="blue"))

        label <- expression(sqrt(U-Z))
        if (ldec) label <- sprintf("%.3f",sqrt((x2-x4)^2+(y2-y4)^2))
        grid.text(label=label,
	      x=xoff+x4+0.4*(x2-x4)+toff*cos(theta),
		  y=yoff+y4+0.4*(y2-y4)-toff*sin(theta),
		  default.units="native",
		  just="left",
		  rot=90-(180/pi)*theta1,		  
		  gp=gpar(col="purple"))
 
} # end of function tdecomp
 


q      <- dat$q
temp   <- tsetup()
width  <- temp$width
height <- temp$height
xoff   <- temp$xoff
yoff   <- temp$yoff
xlim   <- temp$xlim
ylim   <- temp$ylim
maxhw  <- temp$maxhw
rm(temp)

dens <- dat$Nobs
rel  <- dat$rel
res  <- dat$res
bigN <- nrow(dat$res)
unc  <- dat$unc                                     # uncertainty
L    <- dat$L

qstr <- paste("(",
              sprintf("%.2f",q[1]),
	      ", ",
	      sprintf("%.2f",q[2]),
	      ", ",
	      sprintf("%.2f",q[3]),
	      ")",
	      sep="")

sub1 <- bquote(atop(bold(q) * minute == .(qstr)  ,
                   phantom(paste(threshold == .(thresh),phantom(00),
			 delta * p == 1/.(dat$ncirc)))
	           )
	     )

sub2 <- bquote(atop(phantom(bold(q) * minute == .(qstr) ) ,
                   paste(threshold == .(thresh),phantom(00),
			 delta * p == 1/.(dat$ncirc))
	           )
	     )

main1 <- expression(bold(p)[k] - phantom(group("",bold(bar(o)),"|")*bold(p)[k])) 
main2 <- expression(phantom(bold(p)[k]) - group("",bold(bar(o)),"|")*bold(p)[k]) 

tplot(x=q,
        dimnames=c("","",""),
        L=L,
        cex=0,
	col="transparent",
	grid=FALSE,
	bg="white",
        main=main1,
	col.main="black",
	sub=sub1,
	col.sub="blue"
     )
        
tplot(x=q,
        dimnames=c("","",""),
        L=L,
        cex=0,
	col="transparent",
	bg="white",
	grid=FALSE,
        main=main2,
	col.main="red",
	sub=sub2,
	col.sub="black",
	newpage=FALSE
     )
 
     
    pushViewport(viewport(width  = 0.8, 
                          height = 0.8, 
			  xscale = xlim,
                          yscale = ylim))

size <- (maxhw/11)/5

for (i in 1:bigN)
{
   x1 <- xf(matrix(dat$pbin[i,],nrow=1,ncol=3))
   y1 <- yf(matrix(dat$pbin[i,],nrow=1,ncol=3))
   x2 <- xf(matrix(dat$obar[i,],nrow=1,ncol=3))
   y2 <- yf(matrix(dat$obar[i,],nrow=1,ncol=3))

   if (dens[i] > thresh){
	grid.points(x=xoff+x1,
	            y=yoff+y1,
		    pch=21,   # circle
		    size=unit(size,"native"),
		    default.units="native",
		    gp=gpar(col="black",fill="black"))

	grid.points(x=xoff+x2,
	            y=yoff+y2,
		    pch=21,   # circle
		    size=unit(size,"native"),
		    default.units="native",
		    gp=gpar(col="red",fill="transparent"))	
   
        grid.lines(x=xoff+c(x1,x2),
                   y=yoff+c(y1,y2),
		   default.units="native",
		   gp=gpar(col="red"))
		   
		
   }
}
tdecomp( dat,  ldec=TRUE)


if (lsharp) # if we want a sharpness diagram
{
hw <- maxhw  # size of sharpness plot
pushViewport(viewport(
                      x=xlim[2],
                      y=ylim[2],
		      just=c("right","top"),
		      default.units="native",
		      width=hw/2.4,
		      height=hw/2.4,
		      xscale=xlim,
		      yscale=ylim       
		      ))
tsharp(dat=dat,sub="")		      	      
popViewport()	
}


popViewport()

tplot( q,
       L=L,
       col="blue", 
       pch=4, 
       cex=1.1,
       main=paste(""),
       newpage=FALSE,
       bg="transparent",
       grid=FALSE,
       border="transparent")       
}
