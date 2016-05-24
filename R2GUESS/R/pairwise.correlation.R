pairwise.correlation <- function(Y,label.Y=NULL){
par(mar=c(5,5,5,7))
if(is.null(label.Y)) label.Y <- colnames(Y)
colnames(Y) <- label.Y
distance <- cor(Y)
##  Reverse the columns of the matrix so it will be drawn correctly.
ncol <- ncol(distance)
distance2 <- distance[,ncol:1]

#mypalette <- brewer.pal(10,"PiYG")
mypalette <- c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#E6F5D0","#B8E186","#7FBC41"
               ,"#4D9221","#276419")

break.seq <- seq(-1,1,by=0.2)
image(x=1:dim(Y)[2],y=1:1:dim(Y)[2],((distance2)),zlim=c(-1,1),axes=FALSE,xlab="",ylab="",main="Pairwise correlation between the phenotypes",col=mypalette,breaks=break.seq)


par(las=2)
par(cex.axis=1)
axis(2, at = 1:dim(Y)[2],colnames(distance2))
par(cex.axis=1)
axis(1, at = 1:dim(Y)[2],rownames(distance2))
box()
image.plot(x=1:dim(Y)[2],y=1:1:dim(Y)[2],((distance2)),diag=FALSE,zlim=c(-1,1),nlevel=10,xlab="",ylab="",main="Pairwise correlation between the phenotypes",col=mypalette,horizontal = FALSE,legend.only=TRUE,lab.breaks=break.seq,breaks=break.seq)

for(i in 1:(ncol-1)) {
  text(i,1:(ncol-i),round(distance2[i,1:(ncol-i)],digits=2))
}

res <- list(matcor=distance)
}




