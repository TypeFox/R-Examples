dis.cos.cor<-function(fdata1,fdata2=NULL,as.dis=FALSE){
if (is.null(fdata2)) {
     a<-inprod.fdata(fdata1)
     b=diag(a)
     b1=b2=b
#   b<-norm.fdata(fdata1)
#   print(a/(outer(b[,1],b[,1])))

   }
else {
   a<-inprod.fdata(fdata1,fdata2)
   b1<-drop(norm.fdata(fdata1))
   b2<-drop(norm.fdata(fdata2))
#   print(a/(outer(b1[,1],b2[,1])))
}
   if (as.dis) {
   rho<-as.dist(1-abs(a/sqrt(outer(b1,b2))))} else {
   rho<-a/sqrt(outer(b1,b2))}
return(rho)
}
