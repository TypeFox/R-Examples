npmlt <-
function(formula,
         formula.npo=~1,
         random=~1,
         id,
         k=1,
         eps=0.0001,
         start.int=NULL,
         start.reg=NULL,
         start.mp=NULL,
         start.m=NULL,
         link="clogit",
         EB=FALSE,
         maxit=500,
         na.rm=TRUE,
         tol=0.0001)
{

call<-match.call()

# Response variable
a<-model.frame(formula,na.action=na.pass)
resp<-a[,1]
N<-length(resp)

# Design Matrix (Fixed effects)
xp.ii<-model.matrix(formula,data=a)
xp<-xp.ii[,-1]
# Number of Predictors
T<-dim(as.matrix(xp))[2]

# Design Matrix (Random effects)
zp.ii<-model.matrix(random,data=a)
zp<-zp.ii[,-1]
nrp<-dim(as.matrix(zp))[2]

# Design Matrix (Non proportional odds)
npo.design<-model.matrix(formula.npo,data=a)
npo.mtch<-which(noquote(dimnames(xp.ii)[[2]][-1])%in%noquote(dimnames(npo.design)[[2]][-1]))
nnpo<-length(npo.mtch)
#Non % odds indicator
npoind<-array(0,T)
npoind[npo.mtch]<-1

# If id is not provided, use seq(1,N)
if (missing(id)) id<-seq(1,N)

#Select only the complete cases and redefine resp, xp, zp, and id
Temp<-as.data.frame(cbind(resp,id,xp,zp))
if (na.rm==TRUE) Temp<-na.omit(Temp)
if (sum(is.na(Temp))>0) stop("the functions handle only complete cases")

resp<-Temp[,1]
id<-Temp[,2]
if(T>0){xp<-as.matrix(Temp[,(3:(T+2))])}else{xp<-1}
if(nrp>0){zp<-as.matrix(Temp[,((T+3):(T+2+nrp))])}else{zp<-1}

#Stop if data not in the right format
if (min(resp) <= 0) stop("response variable must take values: 1, 2, 3, ...")

# Number of non redundant categories and number of clusters
q<-max(resp)-min(resp)
m<-length(unique(id))

#Number of fixed effectes columns in X matrix
np<-T+nnpo*(q-1)

# Vector with cumulative unique id's
uid<-array(0,m)
for (j in 1:m) uid[j]<-sum(id==unique(id)[j])
cuid<-c(0,cumsum(uid))

# LINK: cumulative logit by default. Use option 'blogit' for baseline logit
linkchoice<-1
if (link=="blogit") linkchoice<-0

# Standard Errors: based on the observed info matrix
SEchoice<-0

# Empirical bayes for random effects: if EB=TRUE EB estimates of the random effects are provided
EBindicator<-0
if (EB==TRUE) EBindicator<-1

# Starting values: the user has the option to provide starting values for none, some or all of the four
# sets of parameters
if (length(start.int)>0 & (!length(start.int)==q))
   stop("start.int is not of the right length")
if (length(start.reg)>0 & (!length(start.reg)==np))
   stop("start.reg is not of the right length")
if (length(start.mp)>0 & (! (length(start.mp)==k*(1+nrp) | length(start.mp)==(k-1)*(1+nrp))) & k>1)
   stop("start.mp is not of the right length")
if (length(start.m)>0 & (!length(start.m)==k) & k>1)
   stop("start.m is not of the right length")

strvlint<-array(0,q)
strvlreg<-NULL
strvlmp<-array(0,k*(1+nrp))
strvlm<-array(0,k)

pred.star.val<-NULL

if (is.null(start.int) | (is.null(start.reg) & T>0)){
   for (j in 1:q){
      temp.resp<-as.numeric(resp<=j)*linkchoice + as.numeric(resp==j)*(1-linkchoice)
      if (T>0) temp<-glm(temp.resp~xp,family=binomial)$coefficients
      if (T==0) temp<-glm(temp.resp~1,family=binomial)$coefficients
      if (is.null(start.int)) strvlint[j]<-temp[1] # to get starting values for the intercepts we do binomial regression on dichotomized y
      if (T>0) pred.star.val<-c(pred.star.val,temp[2:(T+1)])} # and for the coefficients we take the average unless npo
   if (is.null(start.reg) & T>0)
      for (w in 1:T)
         if (w%in%npo.mtch) strvlreg<-c(strvlreg,pred.star.val[seq(w,T*q,T)])
         else strvlreg<-c(strvlreg,mean(pred.star.val[seq(w,T*q,T)]))
}

if (k==1) strvlmp<-array(0,(1+nrp))
if (k>1) strvlmp<-rep(sqrt(2)*gauss.quad(k,"hermite")$nodes[1:(k-1)],(1+nrp))
# mean zero so need first k-1 nodes
# (but its ok to give all of them as they will be ignored by the C program)

weights<-gauss.quad(k,"hermite")$weights
strvlm<-weights/sum(weights)

if (!is.null(start.int)) strvlint<-start.int
if (!is.null(start.reg)) strvlreg<-start.reg
if (!is.null(start.mp) & k>1) strvlmp<-start.mp
if (length(strvlmp)==k*(1+nrp)) strvlmp<-strvlmp[-seq(k,k*(1+nrp),k)]
if (!is.null(start.m) & k>1) strvlm<-start.m

# Length of the output vector
npar<-2*(q+np+(1+1+nrp)*k)+2*(1+nrp)+nrp*(nrp+1)/2+2+2
#q intercepts, np: % and non % predictors, k masses,
#(1+nrp)*k mass points, SE's for all parameters,
#var(RE's') + SE's
#RE's correlation
#but not the errors of the correlations
#logLikelihood, #number of iterations,
#2 flags: cvm and info matrices eigenvalues

#Call to the C function
out<-.C("npmltd",as.integer(resp),as.integer(q),as.integer(N),as.integer(m),as.integer(cuid),as.integer(k),
as.integer(np),as.double(xp),as.double(eps),as.double(strvlint),as.double(strvlreg),as.double(strvlmp),
as.double(strvlm), as.double(array(0,npar)),as.integer(EBindicator),as.double(array(0,(m*(1+nrp)))),
as.integer(linkchoice),as.integer(SEchoice),as.integer(maxit),as.double(zp),as.integer(nrp),
as.double(array(0,(q*N))),as.double(array(0,((q+1)*N))),as.double(tol),as.double(npoind),as.integer(T),
as.double(array(0,((q+np+(nrp+2)*(k-1))^2))))

# Relevant returns of the call to the C function
output<-array(out[[14]][1:npar])
ifelse(nrp>0,eBayes<-matrix(c(out[[16]][1:((1+nrp)*m)]),ncol=(1+nrp),nrow=m),
eBayes<-c(out[[16]][1:((1+nrp)*m)]))
fitted<-matrix(c(out[[22]][1:(q*N)]),nrow=q)
prob<-matrix(c(out[[23]][1:((q+1)*N)]),nrow=(q+1))
CVmat<-matrix(c(out[[27]][1:((q+np+(nrp+2)*(k-1))*(q+np+(nrp+2)*(k-1)))]),nrow=q+np+(nrp+2)*(k-1),ncol=q+np+(nrp+2)*(k-1))

# Give names to the columns of the matrix with the EB estimates
if (nrp>0) dimnames(eBayes)<-list(" "=c()," "=c(dimnames(zp.ii)[[2]]))

#Arrange output
coefficients<-output[1:(q+np)]
if (nrp>0) mass.points<-matrix(output[(q+np+1):(q+np+(1+nrp)*k)],nrow=k,ncol=(1+nrp))
if (nrp==0) mass.points<-output[(q+np+1):(q+np+(1+nrp)*k)]
masses<-output[(q+np+k*(1+nrp)+1):(q+np+k*(1+nrp)+k)]
VCVRE<-output[(q+np+k*(1+nrp)+k+1):(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2)]
ifelse(nrp>0,vcvremat<-diag(VCVRE[1:(1+nrp)]),vcvremat<-VCVRE)
if (nrp>0){
   step<-1
   for (r in 1:nrp){
      for (c in (r+1):(1+nrp)){
         vcvremat[c,r]<-(vcvremat[r,c]<-VCVRE[1+nrp+step])
         step<-step+1}}
}
combmat<-vcvremat
if (nrp>0 & k>1) combmat[upper.tri(combmat)]<-cov2cor(vcvremat)[upper.tri(combmat)]
coefficientsSE<-output[(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+1):(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np)]
if (nrp>0) mass.pointsSE<-matrix(output[(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+1):(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp))],nrow=k,ncol=(1+nrp))
if (nrp==0) mass.pointsSE<-output[(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+1):(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp))]
massesSE<-output[(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+1):(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+k)]
VRESE<-output[(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+k+1):(q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+k+1+nrp)]
LogL<-output[q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+k+1+nrp+1]
iter<-output[q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+k+1+nrp+2]
flagcvm<-output[q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+k+1+nrp+2+1]
flaginfo<-output[q+np+k*(1+nrp)+k+(1+nrp)*(2+nrp)/2+q+np+k*(1+nrp)+k+1+nrp+2+2]

#Output names
cvnames<-dimnames(zp.ii)[[2]]
regsor.names<-NULL
xnames<-dimnames(xp.ii)[[2]][2:(T+1)]
for (i in 1:T) ifelse(npoind[i]==0,regsor.names<-c(regsor.names,xnames[i]),
                      regsor.names<-c(regsor.names,paste(xnames[i],1:q,sep=' ')))
attr(coefficients,'names')<-c(paste('(Intercept)',1:q,sep=' '),regsor.names)[1:(q+np)]
if (nrp>0) dimnames(mass.points)<-list(" "=c(paste('mass point',1:k,sep=' '))," "=c(dimnames(zp.ii)[[2]]))
if (nrp==0) attr(mass.points,'names')<-c(paste('mass point',1:k,sep=' '))
attr(masses,'names')<-c(paste('mass',1:k,sep=' '))
if (nrp>0) dimnames(vcvremat)<-list(" "=c(dimnames(zp.ii)[[2]]), " "=c(dimnames(zp.ii)[[2]]))
if (nrp==0) attr(vcvremat,'names')<-dimnames(zp.ii)[[2]]
if (nrp>0)dimnames(combmat)<-list(" "=c(dimnames(zp.ii)[[2]]), " "=c(dimnames(zp.ii)[[2]]))
if (nrp==0) attr(combmat,'names')<-dimnames(zp.ii)[[2]]
attr(coefficientsSE,'names')<-c(paste('(Intercept)',1:q,sep=' '),regsor.names)[1:(q+np)]
if (nrp>0) dimnames(mass.pointsSE)<-list(" "=c(paste('mass point',1:k,sep=' '))," "=c(dimnames(zp.ii)[[2]]))
if (nrp==0) attr(mass.pointsSE,'names')<-c(paste('mass point',1:k,sep=' '))
attr(massesSE,'names')<-c(paste('mass',1:k,sep=' '))
attr(VRESE,'names')<-paste('SE.Var',dimnames(zp.ii)[[2]],sep=': ')

cv.names<-NULL
if (k>1){
   cv.names<-c(paste('Random.(Intercept)',1:(k-1),sep=''))
   if (nrp>0) for (i in 1:nrp) cv.names<-c(cv.names,paste('Random',cvnames[i+1],1:(k-1),sep=' '))
   cv.names<-c(cv.names,paste('mass point',1:(k-1),sep=' '))
}

dimnames(CVmat)<-list(" "=c(paste('(Intercept)',1:q,sep=''),regsor.names,cv.names),
                      " "=c(paste('(Intercept)',1:q,sep=''),regsor.names,cv.names))

#List of output
fit<-list(call=call,
         formula=formula,
         formula.npo=formula.npo,
         random=random,
         coefficients=coefficients,
         mass.points=mass.points,
         masses=masses,
         vcvremat=vcvremat,
         var.cor.mat=combmat,
         m2LogL=-2*LogL,
         SE.coefficients=coefficientsSE,
         SE.mass.points=mass.pointsSE,
         SE.masses=massesSE,
         VRESE=VRESE,
         CVmat=CVmat,
         eBayes=eBayes,
         fitted=fitted,
         prob=prob,
         nrp=nrp,
         iter=iter,
         maxit=maxit,
         flagcvm=flagcvm,
         flaginfo=flaginfo)
class(fit)<-'npmreg'
return(fit)}
