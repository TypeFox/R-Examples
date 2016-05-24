llikef <-
function(cov1,cov2,dat,m=5,accumulate=TRUE,gradient=FALSE,scale=1,rescale=1,prior=rep(0,length(cov1)+length(cov2))){
	#GRADIENT NOT DEBUGGED FOR SITUATION WITH NO OR MORE THAN 1 coverates
#progression formula, interval censored progression formula
#survival formula, right censored
#combined<-function(dat,cov1,cov2,m,accumulate=TRUE,gradient=TRUE,scale=1){
#dat is a data frame,cov1 and cov2 are formulas(without dep variable) for progression and mortality
#m is the number of intervals, not that max(times)<m
#parms is a m*(m+1)/2 + m +length(cov1)+length(cov2) vector
#outputs a function of parameters (parms)

#Constants that do not depend on parm
dat[,c('tprog0','tprog1','tdeath')]<-rescale*dat[,c('tprog0','tprog1','tdeath')]
ncov1=length(cov1)
ncov2=length(cov2)
nscov1=m*(m+1)/2+m+1 #place of covariates in paramter vector
nscov2=nscov1+ncov1
ms<-dim(dat)[1]
mb=m*(m+1)/2 # number of b parameters
nparms=mb+m+ncov1+ncov2
#s1=attr(
#extract dat
#dtt=data.frame(dat[,c('tprog0','tprog1','tdeath','cdeath')],c1,c2) #needs to be checked
dtt=data.frame(dat[,c('tprog0','tprog1','tdeath','cdeath',cov1,cov2)])
cov1=5:(4+ncov1)
cov2=(5+ncov1):(4+ncov1+ncov2)
#data for cases censored without progression
na1=is.na(dtt$tprog1)
dt1=dtt[na1,] 
dt2=dtt[!na1,] 
tpr0=dt1$tprog0
m11=length(tpr0)
m1=length(dt2[,1])
enum=1:m;
if (m11>0){
mat1=(matrix(floor(tpr0)+1,m11,m)>=matrix(enum,m11,m,byrow=TRUE)) #mask for i...m<= index of last time before progression
mat2=(pmin(matrix(enum,m11,m,byrow=TRUE),matrix(tpr0,m11,m))-matrix(0:(m-1),m11,m,byrow=TRUE))  #min(tprog1,i)-(i-1)
mat1t2=mat1*mat2
}
gm1=diag(1,m,m)
gm2=cbind(matrix(c(0,rep(1,m-1)),ncol=1,nrow=m),matrix(0,m,m-1))
gm3=matrix(0,m,m)
gm3[1,1]=1
if (m>1){
for (i in 2:m) {
	gm1=cbind(gm1,rbind(matrix(0,nrow=i-1,ncol=m-(i-1)),diag(rep(1,m-i+1),m-i+1,m-i+1))) # m x mb
	gm2=rbind(gm2,cbind(matrix(0,nrow=m-i+1,ncol=i-1),c(0,rep(1,m-i)),matrix(0,ncol=m-i,nrow=m-i+1))) # mb x m	
    gm3=rbind(gm3,cbind(matrix(0,nrow=m-i+1,ncol=i-1),c(1,rep(0,m-i)),matrix(0,ncol=m-i,nrow=m-i+1)))  #mb x m
}
}
m1.mb.nparms=c(m1,mb,nparms)
m1.m.nparms=c(m1,m,nparms)
##matrixes umatrices for m=3;
#gm1 
#1 0 0 0 0 0
#0 1 0 1 0 0
#0 0 1 0 1 1
#gm2
#0 0 0 
#1 0 0
#1 0 0
#0 0 0
#0 1 0
#0 0 0
#gm3
#1 0 0 
#0 0 0
#0 0 0
#0 1 0
#0 0 0
#0 0 1
fd=pmin(floor(dt2$tdeath)+1,rep(m,m1))
gm5=((1:m)-1)*(2*(m+1)-(1:m))/2 # 1 x m list of starting places for b vector
gm6=cbind(rep(1:m1,each=m),rep(gm5,times=m1)+rep(fd,each=m)-rep(0:(m-1),times=m1)) #selection matrix for b-vector

#gm5 fd=1,2
#0 3 5
#gm6 will be used to select correct values of b[i,index(fd)] impossible values will be masked out latter files matrix byrow
#1 0+1-0
#1 3+1-1
#1 5+1-2
#2 0+2-0
#2 3+2-1
#2 5+2-2
##gm2a=array(rep(gm2,each=m1),dim=c(m1,mb,m1))
##utr=rep(upper.tri(matrix(1,m,m),diag=FALSE),each=m1,dim(m1,m,m))
##ge4a=array(rep(diag(0:(m-1),m,m),each=m1),dim=c(m1,m,m))
##ge5a=array(rep(gm3,each=m1),dim=c(m1,mb,m))
if (gradient){
	#CODE IS NOT USED numerical gradients worked better
	if (ncov1>0){
	tmps1=array(dtt[,cov1],dim=c(ms,1,ncov1))
	tmps=tmps1
    for (i in (2:m )) {tmps=abind(tmps,tmps1,along=2)}}
   	else tmps=NULL
	if(ncov2>0){
	tmpb1=array(dt2[,cov2],dim=c(ms,1,ncov2)) #NEED TO CHECK WHAT HAPPENS WITH 2 covariat
	tmpb=tmpb1
	for (i in (2:mb)) {tmpb=abind(tmpb,tmpb1,along=2)}}
    else tmpb=NULL
	dlam0=array(c(array(0,dim=c(ms,m,mb)),rep(diag(rep(1,m)),each=ms),tmps),dim=c(ms,m,nparms)) #ms x M x nparms 
	dlam1=array(dlam0[na1,,],dim=c(m11,m,nparms))
	dlam2=array(dlam0[!na1,,],dim=m1.m.nparms)
    db=	array(c(rep(diag(rep(1,mb)),each=m1),array(0,dim=c(m1,mb,m+ncov1)),tmpb),dim=m1.mb.nparms)
	
   gm7= cbind(rep(1:m1,times=m*nparms),rep(rep(gm5,each=m1),times=nparms)+rep(fd,times=m*nparms)-rep(rep(0:(m-1),each=m1),times=m1*nparms),rep(1:nparms,each=m1*m))
	#selection matrix for derivatives	
	#gm7 with m=2 and nparms=2 and m1=2 ready to fill in order of dimension
#1 0+1-0 1
#2 0+2-0 1
#1 0+1-0 1
#2 0+2-0 1
#1 3+1-1 2
#2 3+2-1 2
#1 3+1-1 2
#2 3+2-1 2
}
#data for patients who progressed.
if(m1>0){
mask1=matrix(1:m,m1,m,byrow=TRUE)<=matrix(fd,m1,m)
ppmin=(pmin(matrix(1:m,m1,m,byrow=TRUE),matrix(dt2$tdeath,m1,m))-
                          matrix(0:(m-1),m1,m,byrow=TRUE))*mask1 #  (pmin(tdeath,i)-(i-1)) for i<=fd  
ppmina=(ppmin%*%gm1) #spreads ppmin over mb to get m1 X mb matrix
ppmin1=pmin(matrix(1:m,m1,m,byrow=TRUE),matrix(dt2$tdeath,m1,m)) #min(i,tdeath)			
pm1=pmin(matrix(1:m,m1,m,byrow=TRUE),matrix(dt2$tprog1,m1,m)) # min(i,tprog1)
pm2=pmax(matrix(0:(m-1),m1,m,byrow=TRUE),matrix(dt2$tprog0,m1,m)) #	max((i-1),tprog-)				

#ftmp1=(matrix(floor(dt2$tprog0)+1,nrow=m1,ncol=m)<=matrix(1:m,nrow=m1,ncol=m,byrow=TRUE))&
		#(matrix(1:m,nrow=m1,ncol=m,byrow=TRUE)<=matrix(floor(dt2$tprog1)+1,nrow=m1,ncol=m)) #  floor(tprog0)+1<=i<=floor(tprog1)+1
 #   (matrix(1:m,nrow=m1,ncol=m,byrow=TRUE)<matrix(floor(dt2$tprog1)+1,nrow=m1,ncol=m)) #change 6-4-2012
ft1=rep(dt2$tprog0<dt2$tprog1,m)#interval censored
ft2=pmin(rep(dt2$tprog1,m),rep(1:m,each=m1))>pmax(rep(dt2$tprog0,m),rep(0:(m-1),each=m1))#ith  interval non-empty
#ft3=rep(floor(dt2$tprog0),m)<=rep(0:(m-1),m1,byrow=TRUE) # tprog0 occures on or before i
ft4=rep(floor(dt2$tprog0),m)==rep(0:(m-1),each=m1) # tprog0 in interval i 
ftmp1=matrix((ft1&ft2)|((!ft1)&(ft4)),m1,m)
                                                   
}
#initialization
pdd<-matrix(0,nrow=1,ncol=ms)
if (gradient) dpdd=matrix(0,ms,nparms)
#this is the starting values



#This is the estimation function, it is output from the program in order to have all the onetime work done only one time
outfun=function(parms){
	#print(parms)
if (ncov1>0) parmc=matrix(parms[nscov1:(nscov1+ncov1-1)],ncov1,1);#coefficients for effect on progression of covariates
parml=matrix(parms[(m*(m+1)/2+1):(m*(m+1)/2+m)],1,m)#coeficients for lambda
#This section is for patients who didn't progress  
if(ncov1>0) lam0a=(as.matrix(dtt[,cov1])%*%parmc)%*%matrix(1,ncol=m,nrow=1)+(matrix(1,nrow=ms,ncol=1)%*%parml)  #lam+covariates for prog
else lam0a=(matrix(1,nrow=ms,ncol=1)%*%parml)

if (m11>0) {
lam0=lam0a[na1,]
pdd[na1]=rowSums(-(exp(lam0)*mat1t2))
if (gradient) dpdd[na1,]=apply(c(-exp(lam0)*mat1t2)*dlam1,c(1,3),sum)
}
if (m1>0)
{
#vectorized code for section 2
lam0=lam0a[!na1,]
if (ncov2>0) {
st=5+ncov1
en=st+ncov2-1
b0=as.matrix(dt2[,st:en])%*%matrix(parms[nscov2:(nscov2+ncov2-1)],ncov2,1) #covariates for death m1 x 1
} else b0=matrix(rep(0,m1),m1,1)
parmk=matrix(parms[1:mb],nrow=1,ncol=mb)
b=matrix(parmk,m1,mb,byrow=TRUE)+matrix(b0,m1,mb) #m1 x mb matrix of b's
e1=dt2$cdeath*(matrix(b[gm6],m1,m,byrow=TRUE))+lam0
e2=-(exp(b)*ppmina)%*%gm2
e3=-exp(lam0)%*%upper.tri(matrix(1,m,m),diag=FALSE)
e4=exp(lam0)%*%diag(0:(m-1),m,m)
e5=-ppmin1*exp(b%*%gm3)
ee=e1+e2+e3+e4+e5
u1=exp(b%*%gm3)
u2=exp(lam0)
u=u1-u2


#c1=ezz&(!eze)
#c2=eze
# ett=ifelse(c1,exp(ee)*(pm1-pm2),
#        ifelse(c2,exp(ee+u*pm1),(exp(ee+u*pm1)-exp(ee+u*pm2))/u))
  
  ezz=(u==matrix(0,m1,m))&(!(pm1==pm2)) #u=0 but pm1!=pm2
  eze=(pm1==pm2)
  nezz=(!ezz)&!eze  #ordinary case
  sezz=(sum(ezz)!=0) # at least one case of u=0
  ett=matrix(0,m1,m)

 if (sezz) {
   	ett[ezz]=exp(ee[ezz])*(pm1-pm2)[ezz]
 }	
  tmp11=exp(ee+u*pm1)[nezz]
  tmp12=exp(ee+u*pm2)[nezz]
  tmp1=tmp11-tmp12
  ett[nezz]=tmp1/u[nezz]
  ett[eze]=exp(ee+u*pm1)[eze]
et=rowSums(ett*(ftmp1))	
	pdd[!na1]=log(et)
if (gradient){
  #Code not used possibly wrong
de1=array(dt2$cdeath,dim=m1.m.nparms)*array(db[gm7],dim=m1.m.nparms)+dlam2
de2=apply(array(-c(exp(b)*ppmina)*db,dim=m1.mb.nparms),3,function (x) x%*%gm2)
de3=apply(array(c(-exp(lam0))*dlam2,dim=m1.m.nparms),3,function (x) x%*%upper.tri(matrix(1,m,m),diag=FALSE))
de4=apply(array(c(exp(lam0))*dlam2,dim=m1.m.nparms),3,function (x) x%*%diag(0:(m-1),m,m))
dv=apply(db,3,function(x) x%*%gm3)
de5=c(e5)*dv
dee=de1+array(de2+de3+de4+de5,dim=m1.m.nparms)
du=c(u1)*c(dv)-c(u2)*c(dlam2)           
dett=array(0,dim=m1.m.nparms)
deez=array(ezz,dim=m1.m.nparms)
dett[deez]=ett[ezz]*array(dee,dim=m1.m.nparms)[deez]
tf<-function(x) array(x,dim=m1.m.nparms)[!deez]
if (sum(!deez)>0){                            #some value of ezz is falset
mtmp1=mtmp3=matrix(0,m1,m)


mtmp1[nezz]<- (-tmp1/u[nezz]^2)
mtmp1[eze]<-0
mtmp1[ezz]<-0


mtmp3[nezz]<-(tmp11*pm1[nezz]-tmp12*pm2[nezz])/u[nezz]
mtmp3[eze]<-pm1[eze]*ett[eze]
mtmp3[ezz]<-0

ts1=c(mtmp1)*tf(du)
ts2=c(ett)*tf(dee)
ts3=c(mtmp3)*tf(du)

dett[!deez]=ts1 +ts2+ts3 #may not work if deez ne 0

dett1=array(array(ftmp1,dim=m1.m.nparms)*dett,dim=m1.m.nparms)
dpdd[!na1,]=apply(dett1,c(1,3),sum)/(matrix(et,m1,1)%*%matrix(1,1,nparms))
}
}

}

if (accumulate){

  value=scale*sum(pdd)-0.5*sum(prior*parms[nscov1:nparms]^2)


  if (gradient) attr(value,"gradient")<-scale*colSums(dpdd)
                    -matrix(c(rep(0,nscov1-1),prior*parms[nscov1:nparms]),1,nparms)
   } else {
    value=scale*pdd
	if (gradient) attr(value,"gradient")<-scale*dpdd
    }
return(value)

}
return(outfun)				
}