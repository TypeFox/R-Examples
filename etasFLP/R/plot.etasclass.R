plot.etasclass <-
function(x,pdf=FALSE,file ="etasplot",  ngrid=201,nclass=10,tfixed=0,flag.3D=FALSE,flag.log=FALSE,...){
#  function(x,pdf=FALSE,file ="etasplot",  ngrid=201,flag.3D=FALSE,flag.log=FALSE,ellipse=FALSE,...){
	if(class(x)!="etasclass")stop("argument must be an etasclass object")
	nsimps=1+(ngrid-1)/nclass
        simpson=(nsimps==trunc(nsimps))
	if(!simpson)stop("argument nclass must divide exactly ngrid-1 in order to use Simpson's rule to compute theoretical frequencies")

 scaleres=TRUE
 ellipse=FALSE
# kern.var=x$kern.var
 kern.var=FALSE
# which=1:4
# maxgraph	=7
# ind.graph	=array(FALSE,maxgraph)
# ind.graph[which]=TRUE
 

#tfixed=0 add in the rd file
cat("Computation of the ETAS space intensity integrated intensity on a grid","\n")
cat("for large catalogs computation can take some minutes. Please wait ...","\n")
#####################################################
# plot time intensity # with points
typegraph=1

file=paste(file,".pdf",sep="")
if (pdf){
pdf(file=file,onefile=TRUE)}
#else
#{
#dev.new()
#}
#####################################################
### computation of back intensity and triggered intensity on a given x-y grid
### back.grid computed with the same kernel used for back.dens on observed  points
## back.grid computed on km coordinates

# chamge the is.backconstant switch
magnitude=x$cat$magn1-x$magn.threshold

print(summary(magnitude))

#if(!x$is.backconstant){
{

xcat.km	=x$cat$long
ycat.km	=x$cat$lat
w	=x$rho.weights
hdef	=x$hdef
n	=length(xcat.km)
params	=x$params
        mu= params[1]
        k0= params[2]
        c= params[3]
        p= params[4]
        a= params[5]
        gamma= params[6]
        d= params[7]
        q= params[8]
	#
	    eps=1/n

	    rangex		=range(xcat.km)
	    rangey		=range(ycat.km)
	    eps.x   =diff(rangex)*eps/2
            eps.y   =diff(rangey)*eps/2
            rangex.km  =c(min(rangex)-eps.x,max(rangex)+eps.x)
            rangey.km  =c(min(rangey)-eps.y,max(rangey)+eps.y)

ranget=diff(range(x$cat$time))

alpha	=0
if (kern.var) alpha=hdef[3]

space.grid	=xy.grid(rangex.km,rangey.km,ngrid)
ngridtot	=length(space.grid[,1])
#kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=hdef,kern.var=kern.var,alpha=alpha)
if(!x$is.backconstant){
kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=hdef)
wmat1		=matrix(kde$wmat,n,4)
back.grid	=ranget*mu*kde$z
}
else
{
kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=c(1,1))
wmat1		=matrix(kde$wmat,n,4)
back.grid	=ranget*mu*kde$z^0
}


#### computation MUST be in km IN THIS VERSION because parameters are evaluated in km
ris=.Fortran("etasfull8tintegrated" ,NAOK=TRUE,
			n=as.integer(n),
			mu=as.double(mu),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			a=as.double(a),g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(xcat.km),y=as.double(ycat.km), t=as.double(x$cat$time),m=as.double(magnitude),l=back.grid,
			ngridtot=as.integer(ngridtot), xgrid=as.double(space.grid[,1]),
			ygrid=as.double(space.grid[,2]),
			tmax=as.double(max(x$cat$time)))
			trig.grid	=ris$l
				
### trig.grid intensity on a x-y grid   by time integration of ETAS intensity function
###########################################################################################
}
tot.grid=back.grid+trig.grid

# maps of triggered intensity for a single day

if(tfixed>0){
ind=x$cat$time<tfixed
ris=.Fortran("etasfull8tfixed" ,NAOK=TRUE,
			n=as.integer(sum(ind)),
			mu=as.double(mu),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			a=as.double(a),g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(xcat.km[ind]),y=as.double(ycat.km[ind]), t=as.double(x$cat$time[ind]),m=as.double(magnitude[ind]),l=back.grid,
			ngridtot=as.integer(ngridtot), xgrid=as.double(space.grid[,1]), 
			ygrid=as.double(space.grid[,2]),tfixed=as.double(tfixed))
		        totfixed.grid	=ris$l+back.grid/ranget

if(flag.log)totfixed.grid=log(totfixed.grid)
}

if(flag.log) {
    back.grid=log(back.grid)
    trig.grid=log(trig.grid)
     tot.grid=log( tot.grid)
}
## change grid to degrees

rangex		=range(x$cat.longlat$long)
rangey		=range(x$cat.longlat$lat)
	    eps.x   =diff(rangex)*eps/2
            eps.y   =diff(rangey)*eps/2
            rangex  =c(min(rangex)-eps.x,max(rangex)+eps.x)
            rangey  =c(min(rangey)-eps.y,max(rangey)+eps.y)

space.grid	=xy.grid(rangex,rangey,ngrid)
ngridtot	=length(space.grid[,1])

			x.grid		=seq(rangex[1],rangex[2],length=ngrid)
			y.grid		=seq(rangey[1],rangey[2],length=ngrid)

if(tfixed>0){
### start triggered intensity plotting for a tfixed ########################################
if(!pdf) dev.new()
ind1=x$cat.longlat$time>=tfixed & x$cat.longlat$time<(tfixed+1)

mapxy=map("worldHires",xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),plot=FALSE)

image.plot(x.grid,y.grid,(matrix(totfixed.grid,ngrid,ngrid))
,col=gray.colors(128, start = 0., end = 1., gamma =2 )
,xlab="x-longitude",ylab="y-latitude"
,main=paste("Triggered intensity and observed points at day ",tfixed)
)

grid(col="grey")

map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green"
)
contour(x.grid,y.grid,(matrix(totfixed.grid,ngrid,ngrid)),col="red",add=TRUE)

points(x$cat.longlat$long[ind1],x$cat.longlat$lat[ind1],cex=sqrt(exp(x$cat.longlat$magn1[ind1]))/8,col=4,pch=19)


### end triggered intensity plotting for a tfixed ########################################
}

if(!pdf) dev.new()


### start triggered intensity plotting

mapxy=map("worldHires",xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),plot=FALSE)

image.plot(x.grid,y.grid,(matrix(trig.grid,ngrid,ngrid))
,col=gray.colors(128, start = 0., end = 1., gamma =2 )
,xlab="x-longitude",ylab="y-latitude"
,main="Triggered Intensity"
)

grid(col="grey")

map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green"
)
contour(x.grid,y.grid,(matrix(trig.grid,ngrid,ngrid)),col="red",add=TRUE)

### end triggered intensity plotting

### start background intensity plotting

box()

if(!pdf) dev.new()

image.plot(x.grid,y.grid,(matrix(back.grid,ngrid,ngrid))
,col=gray.colors(128, start = 0., end = 1., gamma =2 ),
xlab="x-longitude",ylab="y-latitude",
main="Background Intensity"
)
      
grid(col="grey")
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green")

contour(x.grid,y.grid,(matrix(back.grid,ngrid,ngrid)),col="red",add=TRUE)

box()

### start total intensity plotting

if(!pdf) dev.new()

image.plot(x.grid,y.grid,(matrix(tot.grid,ngrid,ngrid)),
col=gray.colors(128, start = 0., end = 1., gamma =2 ),
xlab="x-longitude",ylab="y-latitude",main="Total Intensity"
)
      
grid(col="grey")
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green")
contour(x.grid,y.grid,(matrix(tot.grid,ngrid,ngrid)),col="red",add=TRUE)

box()
### start total intensity plotting with observed points
ts=(x$cat$time-min(x$cat$time))/diff(range(x$cat$time))

if(!pdf) dev.new()

image.plot(x.grid,y.grid,
(matrix(tot.grid,ngrid,ngrid)),
col=gray.colors(128, start = 0., end = 1., gamma =2 ),
xlab="x-longitude",ylab="y-latitude",main="Total Intensity with observed points \n Circles area proportional to magnitude; red: recent, blu:older"
)
      
grid(col="grey")
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green")
points(x$cat.longlat$long,x$cat.longlat$lat,cex=sqrt(exp(x$cat.longlat$magn1))/8,col=rgb(ts,0,1-ts),pch=19)
contour(x.grid,y.grid,(matrix(tot.grid,ngrid,ngrid)),col="yellow",add=TRUE)

box()



#####################################################
if((as.numeric(kern.var)*as.numeric(ellipse))==1){

if(!pdf) dev.new()


wx=wmat1[,2]^2
wy=wmat1[,3]^2
wxy=wmat1[,4]*wmat1[,2]*wmat1[,3]
x=xcat.km
y=ycat.km
 plot(x,y,type="p",pch=22,col=1,main="adaptive kernel")
  for (j in 1:n)
  cat(j)
  # internal ellipses
  lines(ellipse(matrix(c(wx[j],wxy[j],wxy[j],wy[j]),2,2),centre=c(x[j],y[j]),level=c(0.500)),col=2)
  # external ellipses
  for (j in 1:n) lines(ellipse(matrix(c(wx[j],wxy[j],wxy[j],wy[j]),2,2),centre=c(x[j],y[j]),level=c(0.999)),col=3)
  }


etas.l=x$l

if(flag.log) etas.l=log(etas.l)

if(flag.3D){

typegraph=2
plot3d(x$cat.longlat$long,x$cat.longlat$lat,etas.l,type="n",zlab=paste("estimated intensity   ","lambda(x,y)"),
xlab="x-longitude",ylab="y-latitude",main="Estimated intensities in observed points")
lines3d(cbind(mapxy$x,mapxy$y,min(etas.l)),col="red")
plot3d(x$cat.longlat$long,x$cat.longlat$lat,etas.l,add=TRUE, type="h")

}
#####################################################################################
### added november 7th, 2014
# computation of  intensities in observed points integrated with respect to time 
#
ngridtot	=n
back.obs	=ranget*mu*kde2dnew.fortran(xcat.km,ycat.km,xcat.km,ycat.km,w=w,factor.xy=1,h=hdef)$z
	
#### computation MUST be in km IN THIS VERSION because parameters are evaluated in km
ris=.Fortran("etasfull8tintegrated" ,NAOK=TRUE,
			n=as.integer(n),
			mu=as.double(mu),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			a=as.double(a),g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(xcat.km),y=as.double(ycat.km), t=as.double(x$cat$time),m=as.double(x$cat$magn1),l=back.obs,
			ngridtot=as.integer(n), xgrid=as.double(xcat.km), ygrid=as.double(ycat.km),tmax=as.double(max(x$cat$time)))
			trig.obs	=ris$l


tot.obs=back.obs+trig.obs

teo1=matrix(0,nclass,nclass)
teo2=matrix(0,nclass,nclass)
teo3=matrix(0,nclass,nclass)
teo4=matrix(0,nclass,nclass)
emp1=teo1
emp2=teo2
emp3=teo3
emp4=teo4

coeff.integr		=simpson.kD(nsimps,2)

	dx      =	diff(rangex.km)/(ngrid-1)
	dy      =	diff(rangey.km)/(ngrid-1)
coeff.integr	=	as.vector(coeff.integr*dx*dy)

x0=min(rangex.km)
y0=min(rangey.km)

xint=diff(rangex.km)/nclass
yint=diff(rangey.km)/nclass

for (i in 1:nclass){
for (j in 1:nclass){

ivert1=(nsimps-1)*(i-1)+1
ivert2=(nsimps-1)*i+1
jvert1=(nsimps-1)*(j-1)+1
jvert2=(nsimps-1)*j+1
xymat=expand.grid(ivert1:ivert2,jvert1:jvert2)
xyvec=ngrid*(xymat[,2]-1)+xymat[,1] 


teovec1=tot.grid[xyvec]
teovec2=back.grid[xyvec]
teo1[i,j]=sum(coeff.integr*teovec1)
#teo2[i,j]=sum(coeff.integr*sqrt(teovec))
teo2[i,j]=sum(coeff.integr*teovec2)

x1=x0+(i-1)*xint
x2=x1+xint
y1=y0+(j-1)*yint
y2=y1+yint
 

indrid=(xcat.km >= x1)&(xcat.km < x2)&(ycat.km >= y1)&(ycat.km < y2)
teo3[i,j]=sum(coeff.integr*sqrt(teovec1))
teo4[i,j]=sum(coeff.integr)
emp1[i,j]=sum(indrid)
emp2[i,j]=sum(x$rho.weights[indrid])
emp3[i,j]=sum(1/sqrt(tot.obs[indrid]))
emp4[i,j]=0
if(emp1[i,j]>0) emp4[i,j]=sum(1/tot.obs[indrid])

}	
}

#####################################################################################
## plotting space residuals 

xx=range(x.grid)
yy=range(y.grid)
hx=diff(xx)/nclass
hy=diff(yy)/nclass

x.grid		=seq(xx[1]+hx/2,xx[2]-hx/2,length=nclass)
y.grid		=seq(yy[1]+hy/2,yy[2]-hy/2,length=nclass)

alpha=255
cblue=  rgb(255:0, 255:0, 255,  maxColorValue = 255,alpha=alpha)
cred =  rgb(255, 255:0, 255:0,  maxColorValue = 255,alpha=alpha)
bluered=c(cred[256:1],cblue)
smax=3

std=(emp1-teo1)/sqrt(teo1)

std4=std<(-smax)
std[std4]=-smax
std4=std>(smax)
std[std4]=smax

if(!pdf) dev.new()

image.plot(x.grid,y.grid,std,zlim=c(-smax,smax),col=bluered,
main="Standardized differences between \n  theoretical and observed frequency (whole model)")

#image.plot(x.grid,y.grid,std,zlim=c(-smax,smax),col=bluered,xaxp=c(xx[1],xx[2],nclass),yaxp=c(yy[1],yy[2],nclass))
grid(nclass)
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",col="black")

############### for background

std=(emp2-teo2)/sqrt(teo2)

std4=std<(-smax)
std[std4]=-smax
std4=std>(smax)
std[std4]=smax

if(!pdf) dev.new()

image.plot(x.grid,y.grid,std,zlim=c(-smax,smax),col=bluered,
main="Standardized differences between  \n theoretical and observed frequency (background only)")

#image.plot(x.grid,y.grid,std,zlim=c(-smax,smax),col=bluered,xaxp=c(xx[1],xx[2],nclass),yaxp=c(yy[1],yy[2],nclass))
grid(nclass)
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",col="black")

############################

######### four graph on a page 

if(!pdf) dev.new()
par(mfcol=c(3,2))

##################### induced seismicity  ############################

#ris=.Fortran("etasfull8reversed" ,NAOK=TRUE,
#			tflag=as.integer(x$onlytime),
#			n=as.integer(n),
#			mu=as.double(mu),k=as.double(k0),
#			c=as.double(c),p=as.double(p),
#			a=as.double(a),g=as.double(gamma),
#			d=as.double(d),q=as.double(q),
#			x=as.double(xcat.km),y=as.double(ycat.km), t=as.double(x$cat$time),m=as.double(magnitude),
#			ltot=as.double(x$l),l=as.double(xcat.km*0))


#	etasn	=ris$l
#y1=log(etasn+0.00000001)


# page7: graph1

#plot(magnitude,y1,xlab="",ylab="log(num. triggered)", main="log(num. triggered) vs. magn1-m0")
#if(var(magnitude>0)){
#    abline(lm(y1 ~ magnitude), lwd = 2, col = 4)
#    lines(lowess(y1 ~ magnitude), lwd = 2, col = 2)
#}
#
## plotting observed vs. observed theoretical
np=200

# page7: graph2
par(cex.main=0.8,mar=c(6, 4, 5, 2) + 0.2)
plot(teo1,emp1,xlab="theor. freq",ylab="empir. freq",main="Total space intens.")
abline(a=0,b=1)
tvec=seq(0,max(teo1),length.out=np)
yvec=1.96*sqrt(tvec)
lines(tvec,tvec+yvec,col=2)
lines(tvec,tvec-yvec,col=2)

##
## plotting space residuals 

xx=range(x.grid)
yy=range(y.grid)
hx=diff(xx)/nclass
hy=diff(yy)/nclass

# page7: graph3

plot(teo2,emp2,xlab="theor. freq",ylab="empir. freq",main="Backgr. space intens.")
abline(a=0,b=1)
tvec=seq(0,max(teo2),length.out=np)
yvec=1.96*sqrt(tvec)
lines(tvec,tvec+yvec,col=2)
lines(tvec,tvec-yvec,col=2)

# page7: graph4

plot(teo1-teo2,emp1-emp2,,xlab="theor. freq",ylab="empir. freq",main="Triggered space intens.")
abline(a=0,b=1)
tvec=seq(0,max(teo1-teo2),length.out=np)
yvec=1.96*sqrt(tvec)
lines(tvec,tvec+yvec,col=2)
lines(tvec,tvec-yvec,col=2)


### page 7 scaled residuals (second column of graph)

# page8: graph1 

## qqplot
#std=(emp1-teo1)/sqrt(teo1)

#qqnorm(std); qqline(std, col = 2)


## plotting residuals vs. theoretical
## on transformed scales
std=(emp1-teo1)/sqrt(teo1)

# page7: graph1

plot(sqrt(teo1),std,xlab="sqrt theor. freq",ylab="std.zd residuals",main="Total space intens.")
abline(a=-1.96,b=0,col=2)
abline(a=+1.96,b=0,col=2)

##
## plotting space residuals 

## plotting residuals vs. theoretical
## on transformed scales
std=(emp2-teo2)/sqrt(teo2)

# page7: graph3

plot(sqrt(teo2),std,xlab="sqrt theor. freq",ylab="std.zd residuals",main="Background space intens.")
abline(a=-1.96,b=0,col=2)
abline(a=+1.96,b=0,col=2)

## plotting residuals vs. theoretical
## on transformed scales
std=(emp1-emp2-(teo1-teo2))/sqrt(teo1-teo2)

# page8: graph4

plot(sqrt(teo1-teo2),std,xlab="sqrt theor. freq",ylab="std.zd residuals",main="Triggered space intens.")
abline(a=-1.96,b=0,col=2)
abline(a=+1.96,b=0,col=2)

############################## time residuals  #######################################
######## 9-12-2014 #############
## transformed points. intensity integrated with respect to space ####
cat("Computation of the time intensity transform","\n")
lambda=mu

ntemp	=n
iprint	=FALSE
ttot	=x$cat$time

rho.tot	=x$rho.s2
t.trasf	=array(0,ntemp)
t.intens=array(0,ntemp)
t.intens[1]=x$back.dens[1]*lambda
ntheta	=ncol(rho.tot)
for (i in 2:ntemp){
times	=ttot[1:i]
range.t	=diff(range(times))
rho.s2	=rho.tot[1:i,]
tmax	=max(times)
magnitudes=magnitude[1:i]
trig=(c+tmax-times)^(-p)
if(p==1){
			 it=log(c+tmax-times)-log(c)
			}
			else
			{
			 it=((c+tmax-times)^(1-p)-c^(1-p))/(1-p)
			}

if(x$onlytime){
		integral=k0*sum(exp(a*magnitudes)*it)
	}
else

##############################################################################
# starting space integration (polar transformation)
##############################################################################
{

m1	=as.vector(exp((a-gamma)*magnitudes))
m2	=as.vector(exp(gamma*magnitudes))
### approximate polar integration to whole space

### approximate polar integration to whole space by division in ntheta triangles centered in xi,yi
etas	=rowSums(m1*m2*((rho.s2*rho.s2/m2+d)^(1-q)-d^(1-q)))

space	=(pi/((1-q)*ntheta))*etas
integral=sum(k0*it*space)
integral.intens=sum(k0*trig*space)
integraltot		=integral	+lambda*range.t*x$back.integral
integraltot.intens	=integral.intens+lambda*x$back.dens[i]

if(iprint) cat(i," ")
t.trasf[i]	=integraltot
t.intens[i]   	=integraltot.intens

}

##############################################################################
#
# end space integration
#
##############################################################################
}

if(!pdf) dev.new()

#close.screen()
par(mfcol=c(2,1))

#split.screen(rbind(c(0,1,0.4,1),c(0,1,0,0.4)))
#screen(1)
plot(t.trasf,type="l",xlab="i",ylab=expression(tau(t[i])),
main="Transformed times")
abline(a=0,b=1,col=2)
#screen(2)
plot(t.intens,type="h",xlab=expression(t[i]),ylab=expression(lambda(t[i])))
close.screen()



################################ end of time residuals


########## end of graphic output ################
if(pdf) dev.off()

return(list(
x.grid=x.grid,
y.grid=y.grid,
space.grid=space.grid,
back.grid=back.grid,
trig.grid=trig.grid,
tot.grid=tot.grid,
back.obs=back.obs,
trig.obs=trig.obs,
tot.obs=tot.obs,
teo1=teo1,
teo2=teo2,
#teo3=teo3,
#teo4=teo4,
emp1=emp1,
emp2=emp2,
#emp3=emp3,
#emp4=emp4,
dx=dx,
dy=dy,
ranget=ranget,
# etasn=etasn,
magnitude=magnitude,
t.trasf=t.trasf
#,wmat
))
}


