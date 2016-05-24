percentile.ple <-
function(dd,p=0.95,gam=0.95,interp=TRUE)  {
# Find Xp the 100pth quantile (percentile) and 100gam%
# nonparametric confidence limit from PLE of F(x)
# USAGE: 
# ARGUMENTS:
# VALUE: Xq the 100q th percentile 
pe<- plekm(dd,gam)
x <- c(0,pe[,1])
#  
qq<-p
p <- c(0,pe$ple)
j <- 0
for( jj in 1:length(x) ) {j<- j+1
   if( p[jj] > qq)  break}

   if(interp){
    xq<- (x[j] - x[j-1])/(p[j] -p[j-1])
    xq<- x[j-1] + xq*( qq - p[j-1])
    }
    else { xq<- (x[j] + x[j-1] )/2  }  

#  UCL    
    j <- 0
    p <- c(0,pe$lower)
   for( jj in 1:length(x) ) {j<- j+1   
   if( p[jj] > qq || j== (length(x)) )  break }
   
   if(interp){
    xqU<- (x[j] - x[j-1])/(p[j] -p[j-1])
    xqU<- x[j-1] + xqU*( qq - p[j-1] )
    }
    else { if(j== length(x)) {xqU<-NA } else{xqU<- x[j] }  }

#  LCL    
    j <- 0
    p<- c(0, pe$upper)
   for( jj in 1:length(x) ) {j<- j+1   
   if( p[jj] > qq)  break}
   if(interp){
    xqL<- (x[j] - x[j-1])/(p[j] -p[j-1])
    xqL<- x[j-1] + xqL*( qq - p[j-1])
    }
    else { xqL <- x[j-1]     }
    
 out<-list(Xp=xq,LXp=xqL,UXp=xqU,p=qq,gam=gam,interp=interp)
out  
}

