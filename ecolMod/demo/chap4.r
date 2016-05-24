##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 4   ##
## Parameterization         ##
##############################

opar <- par()
par(ask=TRUE)
par(mfrow=c(1,1))
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 4  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 4.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

par(mar=c(5.1,4.1,4.1,2.1))

####################################################
# Fig. 4.1
####################################################

ll <- c(0.,1,10,20,40,80,120,160,300,480,700)
pp <- c(0.,1,3,4,6,8,10,11,10,9,8)

plot(ll,pp,xlab= expression("light, muEinst"~ m^{-2}~s^{-1}),
     ylab="production",pch=15,cex=1.5)

fit<-nls(pp ~pmax*2*(1+b)*(ll/iopt)/
                         ((ll/iopt)^2+2*b*ll/iopt+1),            
     start=c(pmax=max(pp),b=0.005,iopt=ll[which.max(pp)]))

summary(fit)

pars <- as.list(coef(fit))

with(pars,
 curve(pmax*2*(1+b)*(x/iopt)/((x/iopt)^2+2*b*x/iopt+1),
       add=TRUE,lwd=2)   )

title(expression (frac(2*pmax*(1+beta)*I/Iopt,
                 (I/Iopt)^2+2*beta*I/Iopt+1)),cex.main=0.8)
subtitle()
########################################
# Figure 4.2. Literature parameters
########################################
 
par(mfrow=c(2,2))
par (las=1)
grow <- Zoogrowth
scoc <- SCOC
ii <- which(grow[,2]>0)
plot(grow[ii,1],grow[ii,2],log="xy",xlab="zooplankton volume, mum3",ylab="" ,
main="maximal growth rate, /hr",pch=16,cex.main=1)

#axis(2,at=c(1e-4,1e-3,1e-2,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
#axis(1,at=c(1e-2,0.1,1,10,100,1000),labels=c("0.01","0.1","1","10","100","1000"))

ll <- lm(log(grow[ii,2])~ log(grow[ii,1]))
rr <- summary(ll)$r.squared
A  <- exp(coef(ll)[1])
B  <- (coef(ll)[2])
x  <- range(grow[ii,1])
lines(x,A*x^B,lwd=2)
AA <- round(A*100)/100
BB <- round(B*100)/100
expr <- substitute(y==A*x^B,list(A=AA,B=BB))
text(100,.0035,expr,adj=0)
expr2 <- substitute(r^2==rr,list(rr=round(rr*100)/100))
text(100,0.002,expr2, adj=0)
writelabel("A")
plot(scoc[,1],scoc[,2],log="xy",xlab="water depth, m",ylab="" ,
     main="SCOC, mmol O2/m2/d",pch=16,xaxt="n",yaxt="n",cex.main=1)

axis(1,at=c(0.5,5,50,500,5000),labels=c("0.5","5","50","500","5000"))
axis(2,at=c(0.1,1,10,100),labels=c("0.1","1","10","100"))

ll <- lm(log(scoc[,2])~ log(scoc[,1]))
rr <- summary(ll)$r.squared
A  <- exp(coef(ll)[1])
B  <- (coef(ll)[2])
x  <- range(scoc[,1],na.rm=TRUE)
lines(x,A*x^B,lwd=2)
AA <- round(A*100)/100
BB <- round(B*100)/100
expr <- substitute(y==A*x^B,list(A=AA,B=BB))
text(1,.1,expr,adj=0)
expr2 <- substitute(r^2==rr,list(rr=round(rr*100)/100))
text(1,0.04,expr2, adj=0)
writelabel("B",at=0.125)

subtitle()

#####################################################
# Figure 4.3. Linearising
#####################################################

par(mfrow=c(2,2))

# the exponential
a<-0.1
b<-2
curve (a*(x^b),0,10,xlab="x",ylab="y")
xx <- seq(0.,10,by=0.25)
noise <-  pmax(0,a*(xx^b)+rnorm(length(xx),sd=0.5) )
text(8,1.5,expression(y==a*x^b))
points(xx,noise)
writelabel("A")

curve(log(a)+b*x,xlim=c(0,2.5),xlab="log(x)",ylab="log(y)")
points(log(xx),log(noise))
writelabel("B")

# The Monod equation
ks <- 1
xx <- seq(0.25,10,by=0.25)
curve (x/(x+ks),0,10,ylim=c(0,1),xlab="x",ylab="y")
noise <-  pmax(0,xx/(xx+ks)+rnorm(length(xx),sd=0.02) )
text(8,0.2,expression(y==frac(x,x+ks)))
points(xx,noise)
writelabel("C")

curve(ks*x+1,xlim=c(0,5),xlab="1/x",ylab="1/y")
points(1/xx,1/noise)
writelabel("D")
subtitle()
#####################################################
# Figure 4.4. Cost surface
#####################################################

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
col=grey(seq(0,0.9,length.out=100))
gg <- outer(seq(3 ,9.5,0.02),seq(-4.8,-1,0.02), FUN=function(x,y) 4*cos(y)+sin(x)-(sin(x)/sqrt(x)*cos(y)*y^2)^2)
persp(gg,col=drapecol(gg,col),border=NA,theta=30,phi=30,axes=TRUE,box=TRUE,
lty=2,xlab="parameter 1",ylab="parameter 2",zlab="cost",ticktype="simple",cex=1.5)
text(0.15,0.22,"2",cex=2)
text(-0.10,0.27,"1",cex=2)
text(-0.15,-0.25,"global minimum")
text(0.1,-0.18,"local minimum")
subtitle()


###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################

########################################
# Figure 4.5. Nonlinear and linear fitting of bioturbation coefficient 
########################################
par(mar=c(5.1,4.1,4.1,2.1))
# the data: 
x <- 0.5:9.5
y <- c(3.9,1.7,1.1,0.5,0.3,0.2,0.1,0.05,0.03,0.02)
plot(y,x,pch=16,ylab="depth, cm",xlab="dpm/cm3",ylim=c(10,0),main="Pb210")

lam <- 0.031

# linear fit
LL  <- lm(log(y)~x)
C0  <- exp(coef(LL)[1])
Db  <- lam/(coef(LL)[2])^2

xx  <-seq(0,10,0.1)
lines(C0*exp(-sqrt(lam/Db)*xx),xx)


fit<-nls(y ~C0*exp(-sqrt(lam/Db)*x),start=c(C0=5,Db=0.1))
C02  <- coef(fit)[1]
Db2  <- coef(fit)[2]
lines(C02*exp(-sqrt(lam/Db2)*xx),xx,lty=2)
legend("topleft",lty=c(1,2),c("linear fit", "nonlinear fit"))


par(new=TRUE,fig=c(0.4,1,0.1,0.7))
plot(y,x,pch=16,ylab="",xlab="",ylim=c(10,0),log="x")
lines(C0*exp(-sqrt(lam/Db)*xx),xx)
lines(C02*exp(-sqrt(lam/Db2)*xx),xx,lty=2)
subtitle()

#####################################################
# Figure 4.6. Price algorithm                        
#####################################################

par(fig=c(0,1,0,1))
gg <- matrix (nr=30,nc=2,data=runif(60))
par(mar=c(1,1,1,1))
plot(gg,xaxt="n",yaxt="n",pch="+",xlim=c(0,1),ylim=c(0,1))
# random select 3 points
#sel <- gg[sample (1:30,3),]  
#points(sel,pch=18,cex=2,col="darkgrey")
sel <- matrix(nr=3,nc=2,data=c(
0.1399560, 0.8592602,0.5327699,
0.52967525,0.77569494, 0.03614811)
)
mirror <- c(0.1665044, 0.8034563)

segments(sel[1,1],sel[1,2],sel[2,1],sel[2,2]) 
segments(sel[3,1],sel[3,2],sel[2,1],sel[2,2]) 
segments(sel[1,1],sel[1,2],sel[3,1],sel[3,2]) 
points(sel,pch=18,cex=2,col="darkgrey")
centroid <- colMeans(sel)
points (centroid[1],centroid[2],pch=19,cex=3,col="darkgrey")
segments(sel[1,1],sel[1,2],centroid[1],centroid[2],lty=2) 
segments(sel[2,1],sel[2,2],centroid[1],centroid[2],lty=2) 
segments(sel[3,1],sel[3,2],centroid[1],centroid[2],lty=2) 

#mirror <- gg[sample(1:30,1),]
points (mirror[1],mirror[2],pch=21,cex=3,col="darkgrey")
newp   <- centroid*2-mirror
points (newp[1],newp[2],pch=17,cex=3,col="black")
arrows(mirror[1],mirror[2],newp[1],newp[2],lty=2,length=0.1)

legend("topright",pch=c(18,19,21,17),c("selected pars","centroid","mirror par", "new par"),
col = c("darkgrey","darkgrey","darkgrey","black"),pt.cex=c(1,2,2,2))
subtitle()

par(mar=c(5.1,4.1,4.1,2.1))
# an application...
amp    <- 6
period <- 5
phase  <- 0.5

x <- sort(runif(20)*13)
y <- amp*sin(2*pi*x/period+phase) +rnorm(20,0,0.05)
plot(x,y,pch=16)

cost <- function(par) sum((par[1]*sin(2*pi*x/par[2]+par[3])-y)^2)
cost(c(1,1,1))
cost(c(6,5,0.5))
p1<-optim(par=c(amplitude=1,phase=1,period=1), cost)
p2<-optim(par=c(amplitude=1,phase=1,period=1), cost,method="SANN")
pp <- optim (par=c(amplitude=1,phase=1,period=1), cost,method="CG")
#            lower=c(0,1e-8,0),upper=c(100,2*pi,100))
p3 <- pricefit(par=c(amplitude=1,phase=1,period=1),minpar=c(0,1e-8,0),
               maxpar=c(100,2*pi,100), func=cost,numiter=3000)

curve(p1$par[1]*sin(2*pi*x/p1$par[2]+p1$par[3]),lty=2,add=TRUE)
curve(p2$par[1]*sin(2*pi*x/p2$par[2]+p2$par[3]),lty=3,add=TRUE)
curve(p3$par[1]*sin(2*pi*x/p3$par[2]+p3$par[3]),lty=1,add=TRUE)

legend ("bottomright",lty=c(1,2,3),c("Price","Mathematical","Simulated annealing"))

subtitle()

#####################################################
# Figure 4.7. Model fit
#####################################################

pgr<-gray.colors(n=25, start=0.95, end=0.0)
with (deepCmin,
filled.contour(x=multser,y=kseries,z=outcost,
             ylab="k (/day)",xlab="multiplication factor (-)",
             main="Model cost landscape",col=pgr,nlevels=25,
             plot.axes={
             axis(1); axis(2);
             points(optpar20$poppar[,2],optpar20$poppar[,1],pch="o",cex=.5);
             points(optpar25$poppar[,2],optpar25$poppar[,1],pch="+",cex=1);
             points(optpar$par[2],optpar$par[1],pch=16,cex=2)
                       }
               )
     )               
subtitle()
par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

