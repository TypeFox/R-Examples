PopGplot <- function(values, colors=FALSE, span=0.1, ylab="", xlab="", ylim=c(min(values,na.rm=TRUE),max(values, na.rm=TRUE))){


if(colors[1]==FALSE){
colors <- rainbow(dim(values)[2])
}

ids <- 1:length(values[,1])
#
# plot first pop
firstpop <- loess(values[,1]~ids, span=span)
plot(predict(firstpop),type="l",xlab=xlab,ylab=ylab,col=colors[1],ylim=ylim,xaxt="n")

if(dim(values)[2]>1){
for(xx in 2:dim(values)[2]){
 popvalues <- values[,xx]
 pop <- loess(popvalues~ids, span=span)
 lines(predict(pop),col=colors[xx])
}

}

}
