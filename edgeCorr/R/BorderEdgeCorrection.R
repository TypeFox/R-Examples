
bec = function(n,pointpattern,Length,Breadth,r){

# applying border edge correction
carrier = array(0,c(n,1))
# right border    
for(i in 1:n){
  carrier[i] = sqrt((pointpattern$x[i] - (Breadth))^{2}  + (pointpattern$y[i] - pointpattern$y[i])^{2})
}
filterright = which(carrier<r)
# left border    
carrierleft = array(0,c(n,1))
for(i in 1:length(carrierleft)){
  carrierleft[i] = sqrt((pointpattern$x[i] - 0)^{2}  + (pointpattern$y[i] - pointpattern$y[i])^{2})
}
filterleft = which(carrierleft<r)
# top border    
carriertop = array(0,c(n,1))
for(i in 1:length(carriertop)){
  carriertop[i] = sqrt((pointpattern$x[i] - pointpattern$x[i])^{2}  + (pointpattern$y[i] - (Length))^{2})
}
filtertop = which(carriertop<r)
# bottom border    
carrierbottom = array(0,c(n,1))
for(i in 1:length(carrierbottom)){
  carrierbottom[i] = sqrt((pointpattern$x[i] - pointpattern$x[i])^{2}  + (pointpattern$y[i] - 0)^{2})
}
filterbottom = which(carrierbottom<r)

combinedfilter = c(filterright,filterleft,filtertop,filterbottom)
globalfilter = unique(combinedfilter)

straysba = cbind(pointpattern$x[globalfilter],pointpattern$y[globalfilter])
correctedba = cbind(pointpattern$x[-globalfilter],pointpattern$y[-globalfilter])

plot(correctedba,col="red",pch=16,cex=1,xlab="x",ylab="y",xlim=c(0,Breadth),ylim=c(0,Length))
points(straysba,col="red",pch=16,cex=2,xlim=c(0,Breadth),ylim=c(0,Length))
abline(h=r,lty=2,col="grey")
abline(v=r,lty=2,col="grey")
abline(v=(Breadth-r),lty=2,col="grey")
abline(h=(Length-r),lty=2,col="grey")

list(globalfilter)

}
