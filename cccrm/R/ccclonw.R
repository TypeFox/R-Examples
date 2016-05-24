ccclonw <-
function(dataset,ry,rind,rtime,rmet,vecD,covar=NULL,rho=0,cl=0.95){

if (length(vecD) == 0) {
stop("Warning: A vector of weights should be provided")}

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


if (length(vecD) != length(unique(dades$time))){
stop("Length of the weight vector must be the number of times")}
D<-diag(vecD)

if ((rho!=0) & (rho!=1)){
stop("Rho must be 0(compound simmetry) or 1 (AR1)")
}

if (rho==0){ 

#Coumpund simmetry model

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
C[,k]<-c(1,contrasts(dades$met)[i,],contrasts(dades$time)[j,],contrasts(dades$met)[i,]*contrasts(dades$time)[j,])
}
}

}


#Weigths matrix


if (nd ==1) auxD<-D

if (nd > 1) {
auxD<-bdiag(list(D,D))
cont<-2
while(cont<nd){
cont<-cont+1
auxD=bdiag(list(auxD,D))
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


alpha=1-cl
Sb<-model.lme$varFix# Var-cov of fixed effects


difmed<-t(L)%*%b
A<-L%*%auxD%*%t(L)

aux1<-(t(difmed)%*%auxD%*%difmed)-sum(diag((A%*%Sb)))
sb<-max(aux1/(nm*(nm-1)),0)
sumd<-sum(D);


# calculating the CCC;
den<-(sumd*(sa+sab+sag+se))+sb
ccc<-(sumd*(sa+sag))/den


# Variance of between-observers variability;

var.sb<-((2*sum(diag(((A%*%Sb)**2))))+(4*(t(b)%*%A%*%Sb%*%A%*%b)))/((nm*(nm-1))^2)



#Covariance between sb and the remeaning parameters;

# dev: Vector of derivatives;

dev.sa<-sumd*(1-ccc)/den
dev.sag<-sumd*(1-ccc)/den
dev.sb<-(-1)*ccc/den
if (sb==0) dev.sb<-0
dev.sab<-sumd*(-1)*ccc/den
dev.se<-sumd*(-1)*ccc/den
dev<-array(c(dev.sag,dev.sab,dev.sa,dev.se,dev.sb),dim=c(1,5))

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

