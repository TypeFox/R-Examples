cccvc2 <-
function(dataset,ry,rind,rmet,covar=NULL,cl=0.95){

options(contrasts=c("contr.treatment","contr.poly"))

dades<-data.frame(dataset)
dades<-rename.vars(dades,from=c(ry,rind,rmet),to=c("y","ind","met"),info=FALSE)
dades$ind<-as.factor(dades$ind)
dades$met<-as.factor(dades$met)
dades$y<-as.numeric(dades$y)


form=y~met
if (length(covar)>0){
form<-as.formula(paste("y~met",paste(covar,sep="+"),sep="+"))}


model.lme<-lme(form,dades,random=list(ind=pdBlocked(list(~1,pdIdent(form=~-1+met)))),method="REML",na.action=na.omit)
if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}
model<-summary(model.lme)

n<-length(unique(dades$ind))# Number fo subjects
k<-length(unique(dades$met))# Number of observers
m<-length(dades$y)/(n*k)

vc<-exp(2*attr(model.lme$apVar,'Pars'))
sa<-vc[1]
sab<-vc[2]
se<-model.lme$sigma^2

S<-4*(vc%*%t(vc))*model.lme$apVar  # Var-cov of variance components

# Building L matrix

aux1<-diag(k)*0
aux1[1,1:k]<-1

for (i in 1:(k-1)) aux1[1+i,1+i]=1

L=array(0,dim=c(k,k*(k-1)/2))
cont=0
for (i in 1:k){
for (j in 1:(k-1)){
if (i>j){
cont=cont+1
 L[,cont]=aux1[,i]-aux1[,j]
}
}
}

# Calculating observers sum of squares
b<-model.lme$coef$fixed[1:k]
difmed=t(L)%*%b
varB<-model.lme$varFix[1:k,1:k]
vardifmed=t(L)%*%varB%*%L
A<-L%*%t(L)
aux2<-(t(difmed)%*%difmed)-sum(diag(A%*%varB))
sb<-max(aux2/(k*(k-1)),0)

# Calculating CCC
den<-sa+sab+sb+se
ccc<-sa/den

# Variance of between-observers variability;

var.sb<-((2*sum(diag(((A%*%varB)**2))))+(4*(t(b)%*%A%*%varB%*%A%*%b)))/((k*(k-1))^2)

#Covariance between sb and the remeaning parameters;

# dev: Vector of derivatives;
alpha=1-cl


dev.sa<-(1-ccc)/den
dev.sb<-(-1)*ccc/den
if (sb==0) dev.sb<-0
dev.sab<-(-1)*ccc/den
dev.se<-(-1)*ccc/den
dev<-array(c(dev.sa,dev.sab,dev.se,dev.sb),dim=c(1,4))

cov.sasb=(-1/n)*(S[1,2]+S[1,3])
cov.sabsb=(-1/n)*(S[2,2]+S[2,3])
cov.sbse=(-1/n)*(S[3,2]+S[3,3])

S2<-array(,c(4,4))
S2[1:3,1:3]<-S
S2[4,]<-c(cov.sasb,cov.sabsb,cov.sbse,var.sb)
S2[1:3,4]<-c(cov.sasb,cov.sabsb,cov.sbse)

varcomp<-c(sa,sab,sb,se)
names(varcomp)<-c("Subjects","Subjects-Method","Method","Error")
est<-ic.ccc(ccc,dev,S2,alpha)

res<-list(ccc=est,vc=varcomp,sigma=S2,model=model)
class(res)<-"ccc"
res 

}


