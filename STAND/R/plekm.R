plekm <-
function(dd,gam=0.95){
#  DESCRIPTION
#
#   Compute product limit estimate estimate (PLE) of F(x) 
#   and 100gam% CLs for left censored data (non-detects) using


# USAGE: 
# plekm(du,gam,xmd))
# ARGUMENTS: 
# du is an n by 2 matrix or data frame
# du[,1]= nonnegative exposure variable (x) in column 1
# du[,2] = det ( 0 for non-detect 1 for detect)
# gam = confidence level for one-sided CLs 
#      ( e.g gam=0.95 generates two-side 90% CLs for ple) 
# xmd = if TRUE and max(x) is non-detect change to detect 
  
# VALUE:   
# data frame with columns
#      
#      x(j) = value of detects
#      ple(j) is PLE of F(x) = Pr[ X <=x ]
#      stder(j) = Standard error of ple
#      lower(j) = LCB for F(x)
#      upper(j) = UCB for F(x)
#      n(j) = number of detects and non-detects <= x(j)
#      d(j)=  number of detects at x(j)
#  NOTE:
#   R function  survreg()is used to calculate Kapalan-Meier estimate 
#   of S(z)  where z = k1 -x .  This technique of reversing the data
#   to obtain right censored data was first suggested by Nelson(1972)
#   conf.type ="plain" is required in survreg() for correct CLs
#
#require(survival)
nx<- dim(dd)[1]  
dd<- dd[ order(dd[,1],dd[,2]),]
#  If xmnd=TRUE and maximum x is a non-detect change it to a detect
# if(dd[nx,2]==0 )if(xmnd){ dd[nx,2]<-1 ;print("maximium x is non-detect") }

if(dd[nx,2]==0 ){ print("maximium x is non-detect") }
 #        if(xmnd) {dd[nx,2]<-1 ;print("CHANGE IT TO A DETECT")}}

# reverse the data by  subtracting x values in column 1 from k1= max(x)

k1<-  max(dd[,1])
d1<- k1 - dd[,1]
dd<- data.frame(d1 ,dd[,2])
cl1<- 1 - (1 -gam)*2
# use next lines test
 ndx.out <-survfit(Surv(dd[,1],dd[,2])~1,se.fit=TRUE,
          data=dd,conf.type="plain",conf.int=cl1) 

ndsout<- summary( ndx.out )

dd[,1] <-  k1 - dd[,1]
# If minimum x is a non-detect add a row to output from survreg()
xmnd<- min( c(dd[ dd[,2]==0,1],k1) ) ; xmd<- min( c(dd[ dd[,2]==1,1]) )

#                      Is smallest x a nondetect ?
if( xmnd < xmd ){      #  yes smallest x is nondetect
x<- k1 - ndsout$time 
nx<- length(x)
a<- c( x,xmnd )
ple<- c(1,ndsout$surv)  ; surv <- 1 - ple

n<- c(ndsout$n.risk,1)
d<- c(ndsout$n.event,0)
lower<-  c(NA,ndsout$lower)
upper<- c(NA,ndsout$upper)
stder<- c(0.0,ndsout$std.err)
out<- cbind(a,ple,stder,lower,upper,n,d)
#          
       }   #  no  smallest x is not a nondetect
else {
a<- k1 - ndsout$time 
nx<- length(a)
n<-ndsout$n.risk
d<-ndsout$n.event

stder<- c( 0.0, ndsout$std.err[1:nx-1] )
ple<- c(1,ndsout$surv[1:nx-1]) 
lower<- c(NA, ndsout$lower[1:nx-1] )
upper<- c( NA, ndsout$upper[1:nx-1] )
out<- cbind( a,ple,stder,lower,upper,n,d)
}
out<- data.frame( out[ order(out[,1]),] )

out
}

