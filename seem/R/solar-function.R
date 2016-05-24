# solar functions

sun.dec <- function(nday){
 dec <- 23.45*sin((2*pi/365)*(nday-81))
 return(dec)
}

sun.elev.max <- function(nday,lat){
 # declination
 dec <- sun.dec(nday)
 if (lat <0) elev <- -dec + 90 - lat
 else elev <- dec + 90 - lat
 return(elev)
}

sun.rad.yr <- function(nday,lat,Lm){
 # declination
 dec <- sun.dec(nday)
 if (lat <0) rad <- Lm - dec*Lm/(90-lat)
 else rad <- Lm + dec*Lm/(90-lat)
 return(rad)
} # end of function

sun.elev.hr <- function(nday,lat,hr.noon){
 # declination
 dec <- sun.dec(nday)
 # hr angle and trig
 nh <- length(hr.noon) 
 hr.angle <- 15*hr.noon
 sin.H <- sin((pi/180)*hr.angle)
 cos.H <- cos((pi/180)*hr.angle)

 # elevation and angle
 sin.elev <- cos((pi/180)*lat)*cos((pi/180)*dec)*cos.H + 
             sin((pi/180)*lat)*sin((pi/180)*dec)
 elev <- asin(sin.elev)*180/pi
 for(i in 1:nh) if(elev[i]<0) elev[i]<-0
 return(elev)
} # end of function

sun.rad.hr <- function(nday,lat,hr.noon,Lm,sdr=0){
 # declination
 elev <- sun.elev.hr(nday,lat,hr.noon)
 rad <- elev
 for(i in 1:length(elev))
  if(elev[i] > 0) rad[i] <- elev[i]*Lm/(90-lat)+ rnorm(1,0,sdr)
 for(i in 1:length(rad))
  if(rad[i] < 0)  rad[i] <- 0 
 return(rad)
} # end of function

# -------------------- hourly rad for several days
sun.rad.hr.mult <- function(nday,lat,Lm,sdr=0,sw.plot=T){

nd <- length(nday)
hr.noon <- seq(-12,+12,0.1); nh <- length(hr.noon)
nl <- max(length(lat),length(Lm))
rad <- matrix(nrow=((nh*nd)-(nd-1)),ncol=nl)

for(i in 1:nl){
   j1<- 1; j2 <- nh 
   for(j in 1:nd){
    rad[j1:j2,i] <- sun.rad.hr(nday[j],lat[i],hr.noon,Lm,sdr)
    j1 <- j2; j2 <- j1+nh-1  
 }
}

hr.cum <- seq(0,24*nd,0.1)
out <- data.frame(hr.cum,rad)

if(sw.plot==T){
 matplot(hr.cum,rad,type="l", col=1, xlab="Hr cumulative",ylab="Radiation [W/m2]")
 abline(h=0,col="grey")
 days <- paste(nday[1])
 for(i in 2:nd) days <- paste(days,nday[i],sep=",")
 mtext(side=3,line=-1,paste("Lat=",lat," Lm=",Lm, " ndays=",days),cex=0.7)
}
return(out)

} # end of function



sun.path <- function(nday,lat){
 # declination
 dec <- sun.dec(nday)
 # hr angle and trig
 hr.noon <- seq(-12,+12,0.1); nh <- length(hr.noon) 
 hr.angle <- 15*hr.noon
 sin.H <- sin((pi/180)*hr.angle)
 cos.H <- cos((pi/180)*hr.angle)

 # elevation and angle
 sin.elev <- cos((pi/180)*lat)*cos((pi/180)*dec)*cos.H + 
             sin((pi/180)*lat)*sin((pi/180)*dec)
 elev <- asin(sin.elev)*180/pi
 # azimuth
 sin.azi <- cos((pi/180)*dec) * sin.H/cos((pi/180)*elev)
 
 # value for testing azimuth
 test.tan <- tan(dec*pi/180)/tan(lat*pi/180)

 kn <- which(hr.noon==0)
 nn <- length(cos.H)
 azi <- array()
 for(k in 1:nn){
  if(cos.H[k] >= test.tan){
   azi[k] <- asin(sin.azi[k])*180/pi
   }
   else{
   if(k <= kn) {corr = -180} else {corr = 180} 
   azi[k] <- corr  - asin(sin.azi[k])*180/pi
   }
}
  for(i in 1:nh) if(elev[i]<0) elev[i]<-0

  return(list(hr.noon=hr.noon,elev=elev, azi=azi))
} # end of function

sun.atmos <- function(nday,lat,rho){
 # Solar constant
 SC = 1.361 # kW/m2
 # ET solar radition I0 kW/m2
 I0 <- SC*(1+0.034*cos((nday)*2*pi/365)) 
 # atmospheric effect
 A <- 1.160 + 0.075 * sin((nday-274)*2*pi/365)
 opt.depth <- 0.174 + 0.035 * sin((nday-100)*2*pi/365)
 elev <- sun.elev.max(nday,lat)
 air.mass <- 1/sin(elev*2*pi/360)
 # Direct
 IB <- I0*exp(-opt.depth*air.mass)
 # diffuse
 IDH <- IB*(0.095 + 0.04*sin((nday-100)*2*pi/365))
 ID <- IDH*(1+cos(pi-elev*2*pi/360))/2
 # reflected
 IBH <- IB*sin(elev*2*pi/360)
 IR <-  rho*(IBH+IDH)*(1+cos(pi-elev*2*pi/360))/2
 # total
 IT <- IB+ID+IR
 I <- cbind(IB,ID,IR,IT)
 return(list(I0=I0,A=A,opt.depth=opt.depth, air.mass=air.mass,I=I))
}


