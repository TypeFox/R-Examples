CVfit <-
function(formula,id,data,family,scale.fix, scale.value,fold,lambda.vec,pindex,eps,maxiter,tol) {

matchedCall <- match.call()
matchedCall[[1]] <- as.name("CVfit")

mf<-model.frame(formula,data)
y<-model.response(mf,"numeric")
X<-model.matrix(formula,data)
#mf$family<-family
#m$link <- family$link
#m$varfun<-family$varfun

if (is.character(family)) family <- get(family)
if (is.function(family))  family <- family()

if(missing(pindex)) pindex=NULL
if(missing(scale.fix))   scale.fix=FALSE
if(missing(scale.value)) scale.value=1
if(missing(eps)) eps=10^-6
if(missing(maxiter)) maxiter=30
if(missing(tol)) tol=10^-3

N<-length(unique(id))
K<-dim(X)[2]-1
nt<-dim(X)[1]/N
nt<-rep(nt,N)

#pay attention to this part  
lam.min<--1
cv.min<-Inf
cv.vect<-NULL
         
for (j in 1:length(lambda.vec))   {
#get one of the lambda's
lam.temp<-lambda.vec[j]
#initial value for cv.values is 0.
cv.value<-0

for (k in 1:fold) {

#select the index that will be omitted.
index.cv<-((k-1)*nt[1]*(N/fold)+1):(k*nt[1]*(N/fold))

#training part#
y.train<-y[-index.cv]
x.train<-X[-index.cv,]
id.train<-id[-index.cv]
       
#compute beta.train
data.train=data.frame("id"=id.train,"y"=y.train,x.train)
mm<-match.call(expand.dots = FALSE)
mm$formula<-formula
mm$id<-id.train
mm$data<-data.train
mm$na.action<-"na.omit"
mm$family<-family
mm$corstr<-"independence"
mm$Mv<-NULL
mm$beta_int<-rep(0,(K+1))  ##NULL##
mm$R<-NULL
mm$scale.fix<-scale.fix
mm$scale.value<-scale.value
mm$lambda<-lam.temp
mm$pindex<-pindex
mm$eps<-eps
mm$maxiter<-maxiter
mm$tol<-tol
mm$silent<-FALSE
mm$lambda.vec<-NULL
mm$fold<-NULL
#mm$link <- NULL
#mm$varfun<NULL

mm[[1]]<-as.name("PGee")

beta.train <- eval(mm, parent.frame())$coefficients

#testing part##
y.cv<-y[index.cv]
x.cv<-X[index.cv,]
id.cv<-id[index.cv]
yy=y.cv
eta=x.cv%*%beta.train 
mu=family$linkinv(eta)
##family$dev.resids gives the square of the residuals
cv.value<-cv.value+sum((family$dev.resids(yy,mu,wt=1)))
} #k
          
cv.vect<-c(cv.vect, cv.value)
           
if(cv.value<cv.min) {
lam.min<-lam.temp
cv.min<-cv.value
}
      
} #j
  
out<-list()
attr(out, "class") <- c("CVfit")
out$fold=fold
out$lam.vect=lambda.vec
out$cv.vect=cv.vect
out$lam.opt=lam.min
out$cv.min=cv.min
out$call <- matchedCall
out   
}
