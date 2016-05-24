cccvc1<-
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


model.lme<-lme(form,dades,random=~1|ind,method="REML",na.action=na.omit)
if(is.character(model.lme$apVar)==TRUE){
stop("Non-positive definite approximate variance-covariance")}
model<-summary(model.lme)





n<-length(unique(dades$ind))# Number of subjects
k<-length(unique(dades$met))# Number of observers
m<-length(dades$y)/(n*k)




# Variance components
tau<-(unlist(model.lme$modelStruct$reStruct))+log(model.lme$sigma)
tau.res<-log(model.lme$sigma)
sa<-exp(2*tau)
se<-(model.lme$sigma^2)

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
den<-sa+sb+se
ccc<-sa/den
names(ccc)<-"CCC"

# Standard error
alpha=1-cl
s.tau<-(model.lme$apVar)[1]# Individual's tau standard error (variance) 
s.res<-(model.lme$apVar)[4]# Error's tau standard error (variance)
s.tau.res<-(model.lme$apVar)[2]#Covariance between  individual and error tau
if ((is.numeric(s.tau)==FALSE) | (is.numeric(s.res)==FALSE) | (is.numeric(s.tau.res)==FALSE)){
print("Convergence Problem")
s.tau<-0
s.res<-0
s.tau.res<-0
}
var.sa<-4*exp(4*tau)*s.tau		#Individual standard error (variance)
var.se<-4*exp(4*tau.res)*s.res	#Error standard error (variance)
cov.sa.se<-4*exp(2*(tau+tau.res))*s.tau.res 	#Covariance between individual and error effects

var.sb<-( 2*sum(diag( ( (A%*%varB)**2 ) ) )+ 4*( t(b)%*%A%*%varB%*%A%*%b ) ) /((k*(k-1))**2)

cov.sa.sb<-(-1/n*m)*cov.sa.se
cov.sb.se<-(-1/n*m)*var.se

S2<-array(,c(3,3))
S2[1,]<-c(var.sa,cov.sa.se,cov.sa.sb)
S2[2,]<-c(cov.sa.se,var.se,cov.sb.se)
S2[3,]<-c(cov.sa.sb,cov.sb.se,var.sb)

dev.sa<-(1-ccc)/den
dev.sb<-(-1)*ccc/den
if (sb==0) dev.sb<-0
dev.se=(-1)*ccc/den

dev<-c(dev.sa,dev.se,dev.sb)

est<-ic.ccc(ccc,t(dev),S2,alpha)
varcomp<-round(c(sa,sb,se),digits=4)
names(varcomp)<-c("Subject","Observer","Random Error")

res<-list(ccc=est,vc=varcomp,sigma=S2,model=model)
class(res)<-"ccc"
res 
}
