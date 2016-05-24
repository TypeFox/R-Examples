cover <-
function(type=c("front","back"),content,col,size,loc,time=NULL){
par(mar=c(1,1,1,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,loc,adj=c(0.5,0.5),cex=size,col=col,labels=content)
if(type=="front"){
rect(30,5,70,15,col="green")
text(50,10,adj=c(0.5,0.5),cex=2,labels="Single click here to start",col="darkgreen")
repeat{ck<-locator(n=1);
if((ck$x>30) & (ck$x<70) & (ck$y>5) & (ck$y<15)) break}}
if(type=="back"){Sys.sleep(time)}
}
