tmap <-
function(dat,     # spatial tverify object
          iyr,             # pick out chosen time of forecast
		  circles=FALSE,   # logical to plot circles
		  fac = 10,        # scale factor for circles
		  theta0=0,        # rotate colour scheme by theta0
		  dich="none",     # show colours as they would appear to a dichromat?
		  m=0.7,           # exponent for saturation assignment
		  palette=TRUE,       # logical to plot triangle legend
		  flip=FALSE)      # flip categories B and A
{ 

lonsize <- (dat$lons[2]-dat$lons[1])/2
latsize <- (dat$lats[2]-dat$lats[1])/2

minlon = min(dat$lons)-lonsize
maxlon = max(dat$lons)+lonsize
minlat = min(dat$lats)-latsize
maxlat = max(dat$lats)+latsize

nlon <- length(dat$lons)
nlat <- length(dat$lats)

tit=paste("Ternary precipitation forecast",sep="")
map("world",xlim=c(minlon,maxlon),ylim=c(minlat,maxlat))

cols <- array(NA,dim=c(nlon,nlat))      
rel  <- array(NA,dim=c(nlon,nlat))  
res  <- array(NA,dim=c(nlon,nlat))     
q    <- array(NA,dim=c(nlon,nlat,3))   

for (ilon in 1:nlon){
   for (ilat in 1:nlat){
         rel[ilon,ilat] <- dat$rel[ilon,ilat]
         res[ilon,ilat] <- dat$res[ilon,ilat]
         q[ilon,ilat,]  <-   dat$q[ilon,ilat,]
     
	 cols[ilon,ilat] <- tcolour(
	                    p=matrix(dat$pred[ilon,ilat,iyr,],ncol=3),
	                    q=matrix(q[ilon,ilat,],ncol=3),
			            m=m,
			            flip=flip,
			            theta0=theta0,
			            dich=dich
			            )

 
			    
	 if (!is.na(cols[ilon,ilat]))
	 {
	    if (circles)
	    {
	       x <- dat$lons[ilon]         
	       y <- dat$lats[ilat]   
               rad <- fac*((sqrt(res[ilon,ilat])-sqrt(rel[ilon,ilat]))/sqrt(res[ilon,ilat])) 	       
               points(x,
	              y,
		      pch=21,
		      cex=rad,
		      bg=cols[ilon,ilat],
		      col="black")    # use score to scale circles 
        }
	    else
	    {	    
	       x <- c(dat$lons[ilon]+lonsize,
	              dat$lons[ilon]-lonsize,
		          dat$lons[ilon]-lonsize,
		          dat$lons[ilon]+lonsize)    # define 4 corners
               y <- c(dat$lats[ilat]+latsize,
	                  dat$lats[ilat]+latsize,
		              dat$lats[ilat]-latsize,
		              dat$lats[ilat]-latsize) 
                polygon(x,y,col=cols[ilon,ilat])  # plot colour if forecast exists
	    }
	  }
   } 
}           

map("world",
    xlim=c(minlon,maxlon),
    ylim=c(minlat,maxlat),
    add=TRUE)

title(main=tit)
if (palette) 
{
    pushViewport(viewport(x=0.9,
                      y=0.2,
		      width=unit(0.2,"snpc"),
		      height=unit(0.2,"snpc")))    # create viewport for legend
   tpalette(q=cbind(1,1,1)/3,
        m=m,
        bars=FALSE,
	flip=flip,
	theta0=theta0,
	dich=dich,
	cex=0.2 
       )     # plot legend     
    popViewport()                   # pop all viewports
}
}
