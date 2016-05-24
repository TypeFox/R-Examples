sc_mult <-
function(X,b,pic=0){
 #calculates scalar multiplication of polygonal fuzzy number X with scalar b
 #X has to be in the same format as the output of the translator function
 ok<-checking(X)
 if(ok==1){
  nl<-nrow(X)
  leftX<-X[1:(nl/2),]
  temp<-X[(nl/2+1):nrow(X),]
  rightX<-temp[nrow(temp):1,]
  if(b>=0){
   sc<-X
   sc$x<-b*X$x
   E<-sc
  }
  if(b<0){
   sc<-X
   temp<-b*X$x
   sc$x<-temp[length(temp):1]
   E<-sc
  }
  if(pic==1){
     limx<-c(min(c(X$x,E$x))-0.25,max(c(X$x,E$x))+0.25)
  plot(X,type="l",xlim=limx, xlab=NA, ylab=expression(alpha),lwd=2, cex.main=1)
  titletxt <- substitute("Fuzzy number and their product by " * b * " (in red)" , list(b = as.character(b)))
      title(main=titletxt,cex.main=1)
 lines(E,type="l",lwd=3,col="red")
   }

invisible(E)
 }
}
