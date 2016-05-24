Vdgraph <-
function(des) 
{
#First check if des is a data frame, and if so convert it to a matrix
  if(is.data.frame(des)) des<-as.matrix(des)
# next check to see if des is centered, and if not center it
  chkm<-apply(des, 2, f)
  if (mean(chkm) != 0) des<-scale(des, scale=FALSE)
  rc<-dim(des)
  ndpts<-rc[1]
   if (ndpts>99) {
         cat("The number of runs in the design matrix is too large","\n")
         } else {
         kvar1<-rc[2]
         }
      check=1+2*kvar1+kvar1*(kvar1-1)/2
       if (ndpts<check) {
         stop("The number of design points will not allow estimation of the quadratic model","\n")
         } else {
         kdv1<-ndpts*kvar1
         }
    if (kvar1>7) {
         stop("The maximum number of independent variables allowed is 7","\n")
                }		 
				rdes<-des
    dim(rdes)<-ndpts*kvar1
cat("number of design points=",ndpts,"\n")
cat("number of factors=",kvar1,"\n")
v<-Vardsgr(ndpts,kvar1,kdv1,rdes)
vdgr<-matrix(v,ncol=4)
colnames(vdgr)<-c("Radius","Maximum","Minimum","Average")
# Fix fortran bug
     icnt<-0
	 for (k in 1:21) {
	  if (vdgr[k,3]==vdgr[k,4] & vdgr[k,2]==vdgr[k,4]) {
	       icnt<-icnt+1}
		   else { if (icnt>=4){
		        vdgr[k,2]<-vdgr[k,4]}}
				      }
#modifies plot region
par(mai=c(1.0,.75,1.0,.25))
#Creates the plot region
maxy=max(vdgr[,2])
aminy<-min(vdgr[,3])
miny=0
maxx=sqrt(kvar1)
#cat("maxy=",maxy,"\n")
#cat("aminy=",aminy,"\n")
#cat("miny=",miny,"\n")
#cat("maxx=",maxx,"\n")
plot(c(0,maxx),c(miny,maxy),type="n",ylab="Scaled Variance",xlab="Radius",
     main="Variance Dispersion Graph")
#Adds the line for maximum variance
lines(vdgr[,1],vdgr[,2],lty=2,col="red")
#Adds the line for minimum variance
lines(vdgr[,1],vdgr[,3],lty=4,col="blue")
#Adds the line for average variance
lines(vdgr[,1],vdgr[,4],lty=1,col="black")
#Adds the legend
legend("topleft",inset=.02,legend=,c("Max","Min","Avg"),lty=(c(2,4,1)),col=(c("red","blue","black")))
#gets default graphic parameters
defpar<-par()
return(vdgr)
}

