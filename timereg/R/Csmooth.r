CsmoothB<-function(cumB,xval,b,degree=1,coef=1) 
{
cumB<-as.matrix(cumB); 
if(is.matrix(cumB) == TRUE) px<-as.integer(dim(cumB)[2]);
if(is.matrix(cumB) == TRUE) nx<-as.integer(dim(cumB)[1]);
nval<-length(xval); 
bhat<-matrix(0,nval,px-1);
bhat<-cbind(xval,bhat)
b<-matrix(b,nval,px-1); 
#print(b); 

sout<-.C("smoothB", 
as.double(cumB),as.integer(nx),as.integer(px),
as.double(bhat),as.integer(nval),
as.double(b),as.integer(degree),as.integer(coef),PACKAGE="timereg")

bhat<-matrix(sout[[4]],nval,px);
return(bhat)
}


Csmooth2B<-function(cumB,xval,b,degree=1,coef=1) 
{
cumB<-as.matrix(cumB); 
if(is.matrix(cumB) == TRUE) px<-as.integer(dim(cumB)[2]);
if(is.matrix(cumB) == TRUE) nx<-as.integer(dim(cumB)[1]);
nval<-length(xval); 
bhat<-matrix(0,nval,px-1);
bhat<-cbind(xval,bhat)
b<-matrix(b,nval,px-1); 

sout<-.C("smooth2B", 
as.double(cumB),as.integer(nx),as.integer(px),
as.double(bhat),as.integer(nval),
as.double(b),as.integer(degree),as.integer(coef),PACKAGE="timereg")

bhat<-matrix(sout[[4]],nval,px);
return(bhat)
}

#bhat<-smoothB(xval[100],3,SAtime$cum)
#ud<-CsmoothB(SAtime$cum,xval,3,degree=2)[100,];
#ud<-Csmooth2B(SAtime$cum,xval,bx,degree=2);
#(bhat/ud)
#bx<-c(rep(2,90),rep(3,211))

localTimeReg<-function(times,response,designX,xval,b,lin=1,dens=0) 
{
designX<-as.matrix(designX); 
if(is.matrix(designX) == TRUE) px<-as.integer(dim(designX)[2]);
if(is.matrix(designX) == TRUE) nx<-as.integer(dim(designX)[1]);
nval<-length(xval); 

bhat<-matrix(0,nval,(lin+1)*px);
bhat<-cbind(xval,bhat);
fdens<-cbind(xval,0); 

b<-matrix(b,nval,1); 

sout<-.C("localTimeReg", 
as.double(designX),as.integer(nx),as.integer(px),
as.double(times),as.double(response),
as.double(bhat),as.integer(nval),
as.double(b),as.integer(lin),as.double(fdens),PACKAGE="timereg")

bhat<-matrix(sout[[6]],nval,(lin+1)*px+1);
if (lin==2) bhat[,4]<-bhat[,4]*2
if (lin==3) bhat[,5]<-bhat[,5]*6

fdens<-matrix(sout[[10]],nval,2);
colnames(fdens)<-c("density","Ddensity"); 
if (dens==1) bhat<-cbind(bhat,fdens); 

return(bhat)
}
