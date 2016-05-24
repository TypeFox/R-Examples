
pricefit <- function (
            par,                          # initial par estimates
            minpar=rep(-1e8,length(par)), # minimal parameter values
            maxpar=rep(1e8,length(par)),  # maximal parameter values
            func,                         # function to minimise
            npop=max(5*length(par),50),   # nr elements in population
            numiter=10000,                # number of iterations
            centroid = 3,                 # number of points in centroid 
            varleft  = 1e-8,              # relative variation upon stopping
            ...)

{

# Initialisation

 cost  <- function (par) func(par,...)
 npar  <- length(par)
 tiny  <- 1e-8
 varleft<-max(tiny,varleft)
 
 

 populationpar  <- matrix(nrow=npop,ncol=npar,byrow=TRUE,
             data= minpar+runif(npar*npop)*rep((maxpar-minpar),npop))
 colnames(populationpar)<-names(par)
 populationpar[1,]<-par

 populationcost <- apply(populationpar,FUN=cost,MARGIN=1)
 iworst         <- which.max(populationcost)
 worstcost      <- populationcost[iworst]


# Hybridisation phase   
 iter<-0
 while (iter<numiter & (max(populationcost)-min(populationcost))
                                  >(min(populationcost)*varleft))
 {
   iter<-iter+1
   
   selectpar <- sample(1:npop,size=centroid)  # for cross-fertilisation
   mirrorpar <- sample(1:npop,size=1)         # for mirroring  

   newpar    <- colMeans(populationpar[selectpar,])    # centroid
   newpar    <- 2*newpar - populationpar[mirrorpar,]   # mirroring
     
   newpar    <- pmin( pmax(newpar,minpar) ,maxpar)
                        
   newcost   <- cost(newpar)


    if (newcost < worstcost)  
     { 
       populationcost[iworst] <-newcost
       populationpar [iworst,]<-newpar  
       iworst     <- which.max(populationcost) # new worst member
       worstcost  <- populationcost[iworst]
     }
  } # end j loop


  ibest    <- which.min(populationcost)
  bestpar  <- populationpar[ibest,]
  bestcost <- populationcost[ibest]
return (list(par = bestpar, cost = bestcost, 
             poppar = populationpar, popcost=populationcost))

}


################################################
# drawing of a DILUTION CULTURE
################################################


dilution <- function (main=c("",""),int="")
{
emptyplot()
dx <- 0.01
dx2<- 0.02
dx3<- dx2+dx

xx <- c(0.30,0.8,0.95)
yy <- c(0.1,0.65,0.45)

dxb<-0.1
dx4<-0.5
col <- grey(0.9)

lines(c(xx[2],xx[2],xx[1],xx[1],0.35-dx),c(yy[1]+dx4-dx3,yy[1],yy[1],yy[2],yy[2]))
lines(c(0.35+dx3,xx[2],xx[2]),c(yy[2],yy[2],yy[1]+dx4+dx))
rect(xx[1]+dx,yy[1]+dx,xx[2]-dx,yy[2]-dxb ,lwd=1,col=col,border=NA)

# protrusion
rect(xx[2]-dx,yy[1]+dx4,xx[3],yy[1]+dx4-dx2,col=col,border=NA)
rect(xx[3]-dx2,yy[1]+dx4,xx[3],yy[3],col=col,border=NA)
segments(xx[2],yy[1]+dx4+dx,xx[3]+dx,yy[1]+dx4+dx)
segments(xx[2],yy[1]+dx4-dx3,xx[3]-dx3,yy[1]+dx4-dx3)
segments(xx[3]-dx3,yy[1]+dx4-dx3,xx[3]-dx3,yy[3])
segments(xx[3]+dx,yy[1]+dx4+dx,xx[3]+dx,yy[3])

filledcircle(mid=c(mean(c(xx[3]+dx,xx[3]-dx3)),yy[3]-dx2),r1=0.01,lwd=1,col=col)
filledcircle(mid=c(mean(c(xx[3]+dx,xx[3]-dx3)),yy[3]-2.5*dx2),r1=0.01,lwd=1,col=col)

text(mean(xx[1:2]),yy[2]+dx2,main[2])
p1 <-c(mean(xx[1:2]),mean(yy[1:2]))

arrows(xx[2]+dx,yy[1]+dx4+2*dx, xx[3]+dx,yy[1]+dx4+2*dx,lwd=2,length=0.1)
text(mean(c(xx[2]+dx,xx[3]+dx)),yy[1]+dx4+4*dx,int,font=2)
dxb<- 0.03
dx4<- 2*dx+dx2
x2 <- c(0.05,0.25,0.35)
y2 <- c(0.70,0.94,0.68)
col <- grey(0.8)

lines(c(x2[2],x2[1],x2[1],x2[2],x2[2]),c(y2[1],y2[1],y2[2],y2[2],y2[1]+dx4))
rect(x2[1]+dx,y2[1]+dx,x2[2]-dx,y2[2]-dxb ,lwd=1,col=col,border=NA)
rect(x2[2]-dx,y2[1]+dx,x2[3],y2[1]+dx3,col=col,border=NA)
rect(x2[3],y2[1]+dx3,x2[3]+dx2,y2[3],col=col,border=NA)
# tube
segments(x2[2]-dx,y2[1],x2[3]-dx,y2[1])
segments(x2[2],y2[1]+dx4,x2[3]+dx3,y2[1]+dx4)
segments(x2[3]+dx3,y2[1]+dx4,x2[3]+dx3,y2[3])
segments(x2[3]-dx,y2[1],x2[3]-dx,y2[3])

filledcircle(mid=c(mean(c(x2[3]+dx3,x2[3]-dx)),y2[3]-dx2),r1=0.01,lwd=1,col=col)
filledcircle(mid=c(mean(c(x2[3]+dx3,x2[3]-dx)),y2[3]-2.5*dx2),r1=0.01,lwd=1,col=col)


text(mean(x2[1:2]),y2[2]+dx2,main[1])
p2 <-c(mean(x2[1:2]),mean(y2[1:2]))

arrows(x2[2]+dx,y2[1]+dx4+dx, x2[3]+dx3-dx,y2[1]+dx4+dx,lwd=2,length=0.1)
text(mean(c(x2[2]+dx,x2[3]+dx3-dx)),y2[1]+dx4+3*dx,int,font=2)
list(p1=p1,p2=p2)
}

