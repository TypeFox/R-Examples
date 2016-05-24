# ex2.10.R
require(sde)
set.seed(123)
W <- vector(14,mode="list")
W[[1]] <- BBridge(1,1,0,1,N=2)

for(i in 1:13){
 cat(paste(i,"\n"))
 n <- length(W[[i]])
 t <- time(W[[i]])
 w <- as.numeric(W[[i]]) 
 tmp <- w[1]
 for(j in 1:(n-1)){  
   tmp.BB <- BBridge(w[j],w[j+1],t[j],t[j+1],N=2)
   tmp <- c(tmp, as.numeric(tmp.BB[2:3]))
  }
  W[[i+1]] <- ts(tmp,start=0,deltat=1/(2^(i+1)))
}
min.w <- min(unlist(W))-0.5
max.w <- max(unlist(W))+0.5
opar <- par(no.readonly = TRUE)
par(mfrow=c(7,2),mar=c(3,0,0,0))
for(i in 1:14){
 plot(W[[i]], ylim=c(min.w, max.w),axes=F)
 if(i==1)
  axis(1,c(0,0.5,1))
 if(i==2)
 axis(1,c(0,0.25,0.5,0.75,1))
 if(i>2)		
  axis(1,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
 text(0.5,2.2,sprintf("N = %d",2^i))
}
par(opar)
