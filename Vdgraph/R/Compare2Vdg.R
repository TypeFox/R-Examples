Compare2Vdg <-
function(des,des2,name1=deparse(substitute(des)),name2=deparse(substitute(des2)),ncolleg=1) 
{
#Check to see of ncolleg is 1 or 2
if(!(ncolleg !=1 | ncolleg !=2))
  stop("ncolleg, the number of columns in the legend must be set equal to 1 or 2","\n")
#Check if names are too long for the legend
cn1<-nchar(name1)
cn2<-nchar(name2)
mcn<-max(cn1,cn2)
  if(mcn>40) {
    stop("The name character strings are too long to fit in the legend","\n")
             }
#Next check if des is a data frame, and if so convert it to a matrix
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
       if (ndpts<check){
         stop("The number of design points will not allow estimation of the quadratic model","\n")
         } else {
         kdv1<-ndpts*kvar1
         }
    if (kvar1>7) {
         stop("The maximum number of independent variables allowed is 7","\n")
                }		 
				rdes<-des
    dim(rdes)<-ndpts*kvar1
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
nfac1<-kvar1

#First check if des2 is a data frame, and if so convert it to a matrix
  if(is.data.frame(des2)) des2<-as.matrix(des2)
# next check to see if des2 is centered, and if not center it
  chkm<-apply(des2, 2, f)
  if (mean(chkm) != 0) des2<-scale(des2, scale=FALSE)
  rc<-dim(des2)
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
    if (kvar1 != nfac1) {
         stop("The number of factors in the two designs must match","\n")
                }		 
	 
				rdes<-des2
    dim(rdes)<-ndpts*kvar1
v<-Vardsgr(ndpts,kvar1,kdv1,rdes)
vdgr2<-matrix(v,ncol=4)
colnames(vdgr2)<-c("Radius","Maximum","Minimum","Average")
# Fix fortran bug
     icnt<-0
	 for (k in 1:21) {
	  if (vdgr2[k,3]==vdgr2[k,4] & vdgr2[k,2]==vdgr2[k,4]) {
	       icnt<-icnt+1}
		   else { if (icnt>=4){
		        vdgr2[k,2]<-vdgr2[k,4]}}
				      }

#modifies plot region
par(xpd=TRUE,mai=c(2.5,.75,.5,.25))
#Creates the plot region
maxy1=max(vdgr[,2])
maxy2=max(vdgr2[,2])
aminy1<-min(vdgr[,3])
aminy2<-min(vdgr2[,3])
maxy<-max(maxy1,maxy2)
miny=0
ledy<-maxy/3.8
maxx=sqrt(kvar1)
plot(c(0,maxx),c(miny,maxy),type="n",ylab="Scaled Variance",xlab="Radius",
     main="Variance Dispersion Graph")
#Adds the line for maximum variance of des
lines(vdgr[,1],vdgr[,2],lty=2,col="royalblue")
#Adds the line for minimum variance of des
lines(vdgr[,1],vdgr[,3],lty=4,col="darkblue")
#Adds the line for average variance of des
lines(vdgr[,1],vdgr[,4],lty=1,col="blue")
#Adds the line for maximum variance of des2
lines(vdgr[,1],vdgr2[,2],lty=2,col="red")
#Adds the line for minimum variance of des2
lines(vdgr[,1],vdgr2[,3],lty=4,col="darkred")
#Adds the line for average variance of des2
lines(vdgr[,1],vdgr2[,4],lty=1,col="magenta")
# makes names for legend
Xname1<-paste("Max(",name1,")",sep="")
Nname1<-paste("Min(",name1,")",sep="")
Aname1<-paste("Avg(",name1,")",sep="")
Xname2<-paste("Max(",name2,")",sep="")
Nname2<-paste("Min(",name2,")",sep="")
Aname2<-paste("Avg(",name2,")",sep="")

##Adds the legend
legend(0,-ledy,legend=c(Xname1,Nname1,Aname1,Xname2,Nname2,Aname2),lty=(c(2,4,1,2,4,1)),
col=(c("royalblue","darkblue","blue","red","darkred","magenta")),ncol=ncolleg)

}

