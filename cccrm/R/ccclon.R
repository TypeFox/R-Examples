ccclon <-
function(dataset,ry,rind,rtime,rmet,covar=NULL,rho=0,cl=0.95){

dades<-data.frame(dataset)
dades<-rename.vars(dades,from=c(ry,rind,rmet,rtime),to=c("y","ind","met","time"),info=FALSE)
dades$ind<-as.factor(dades$ind)
dades$met<-as.factor(dades$met)
dades$time2<-as.numeric(dades$time)
dades$time<-as.factor(dades$time)
dades$y<-as.numeric(dades$y)

form=y~met+time+met*time
if (length(covar)>0){
form<-as.formula(paste("y~met+time+met*time",paste(covar,sep="+"),sep="+"))}


if ((rho!=0) & (rho!=1)){
stop("Rho must be 0(compound simmetry) or 1 (AR1)")
}


if (rho==0){ 

#Compound simmetry model

model.lme<-lme(form,dades,random=list(ind=pdBlocked(list(~1,pdIdent(form=~-1+met),pdIdent(form=~-1+time)))))
if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}
model<-summary(model.lme)

# Variance components



vc<-exp(2*attr(model.lme$apVar,'Pars'))
sa<-vc[1]
sab<-vc[2]
sag<-vc[3]
se<-model.lme$sigma^2

S<-4*(vc%*%t(vc))*model.lme$apVar  # Var-cov of variance components
}

if (rho==1){


#AR1 model

model.lme<-lme(form,dades,random=list(ind=pdBlocked(list(~1,pdIdent(form=~-1+met),pdIdent(form=~-1+time)))),correlation=corAR1(form=~time2|ind/met))
if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}
model<-summary(model.lme)

# Variance components

vc<-exp(2*attr(model.lme$apVar,'Pars'))[c(1:3,5)]
sa<-vc[1]
sab<-vc[2]
sag<-vc[3]
se<-model.lme$sigma^2
S<-4*(vc%*%t(vc))*model.lme$apVar[c(1:3,5),c(1:3,5)]  # Var-cov of variance components


}
# Dimensions
ns<-length(unique(dades$ind))
nt<-length(unique(dades$time))
nm<-length(unique(dades$met))


if ((rho==0) | (rho==1)){

b<-as.matrix(model.lme$coef$fixed)

# Design matrix
nd<-nm*(nm-1)/2
C<-array(0,dim=c(length(b),nt*nm))
k<-0
for (i in 1:nm){
for (j in 1:nt){
k<-k+1
C[,k]<-c(1,contrasts(dades$met)[i,],contrasts(dades$time)[j,],c(contrasts(dades$met)[i,]%*%t(contrasts(dades$time)[j,])))
}
}




# Difference between methods matrix
L<-array(0,dim=c(length(b),nt*nd))
k<-0
for (i in 1:(nt*(nm-1))){
for (j in (1:(nm-1))){
if ((i+nt*j)<=(nt*nm)){
k<-k+1
L[,k]=C[,i]-C[,i+nt*j]
}
}
}



Sb<-model.lme$varFix# Var-cov of fixed effects


difmed<-t(L)%*%b
A<-L%*%t(L)

aux1<-(t(difmed)%*%difmed)-sum(diag((A%*%Sb)))
sb<-max(aux1/(nm*(nm-1)*nt),0)
}
# calculating the CCC;
den<-sa+sab+sag+se+sb
ccc<-(sa+sag)/den



# Variance of between-observers variability;

var.sb<-((2*sum(diag(((A%*%Sb)**2))))+(4*(t(b)%*%A%*%Sb%*%A%*%b)))/((nm*(nm-1)*nt)^2)

#Covariance between sb and the remeaning parameters;

# dev: Vector of derivatives;
alpha=1-cl


dev.sa<-(1-ccc)/den
dev.sag<-(1-ccc)/den
dev.sb<-(-1)*ccc/den
if (sb==0) dev.sb<-0
dev.sab<-(-1)*ccc/den
dev.se<-(-1)*ccc/den
dev<-array(c(dev.sa,dev.sab,dev.sag,dev.se,dev.sb),dim=c(1,5))

cov.sasb=(-1/ns)*(S[1,2]+S[1,4])
cov.sabsb=(-1/ns)*(S[2,2]+S[2,4])
cov.sagsb=(-1/ns)*(S[3,2]+S[3,4])
cov.sbse=(-1/ns)*(S[4,2]+S[4,4])

S2<-array(,c(5,5))
S2[1:4,1:4]<-S
S2[5,]<-c(cov.sasb,cov.sabsb,cov.sagsb,cov.sbse,var.sb)
S2[1:4,5]<-c(cov.sasb,cov.sabsb,cov.sagsb,cov.sbse)

varcomp<-c(sa,sab,sag,sb,se)
names(varcomp)<-c("Subjects","Subjects-Method","Subjects-Time","Method","Error")
est<-ic.ccc(ccc,dev,S2,alpha)

res<-list(ccc=est,vc=varcomp,sigma=S2,model=model)
class(res)<-"ccc"
res 

}