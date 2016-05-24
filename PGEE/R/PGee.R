PGee <-
function(formula,id,data,na.action=NULL,family=gaussian(link="identity"),
corstr="independence",Mv=NULL,beta_int=NULL,R=NULL,scale.fix=FALSE,
scale.value=1,lambda,pindex=NULL,eps=10^-6,maxiter=30,tol=10^-3,silent=FALSE)  {

call <- match.call()
m <- match.call(expand.dots = FALSE)

m$beta_int <- m$family <- m$link <- m$varfun<-
m$corstr <- m$Mv<- m$R <-
m$scale.fix <- m$scale.value <-
m$lambda <-m$eps<-m$pindex<-
m$maxiter <- m$tol <-m$silent <- NULL

if(is.null(m$id)) m$id<-as.name("id")

if(!is.null(m$na.action) && m$na.action != "na.omit") {
warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
m$na.action <- as.name("na.omit")
}

m[[1]] <- as.name("model.frame")
m <- eval(m, parent.frame())
Terms <- attr(m, "terms")
y <- model.extract(m, "response")
X<- model.matrix(Terms, m, contrasts)

xnames <- dimnames(X)[[2]]
if(is.null(xnames)) {
xnames <- paste("x", 0:K, sep = "")
dimnames(X) <- list(NULL, xnames)
}

id<-model.extract(m, id)

if(is.null(id)) {
stop("Id variable not found!")
}

if(length(id) != length(y))  stop("Id and y do not have the same length!")

if(!(is.double(X)))  X <- as.double(X)
if(!(is.double(y)))  y <- as.double(y)
if(!(is.double(id))) id <- as.double(id)

N<-length(unique(id))
K<-ncol(X)-1

avec <- as.integer(unlist(lapply(split(id, id), "length")))
maxclsz <-max(avec)
maxcl <- maxclsz
nt<-avec
nobs<-sum(nt)

if(!(is.double(N)))      N <- as.double(N)
if(!(is.double(maxcl)))  maxcl <- as.double(maxcl)
if(!(is.double(nobs)))   nobs <- as.double(nobs)

if(missing(lambda)) stop("A value is not assiged for lambda!")

if(missing(pindex)) pindex=NULL

if(missing(family)) family=gaussian(link="identity")

if(missing(corstr)) corstr="independence"

if(missing(Mv)) Mv<-NULL

if(corstr=="stat_M_dep" && is.null(Mv)) stop("corstr is assumed to be 'stat_M_dep' but Mv is not specified!")

if(corstr=="non_stat_M_dep" && is.null(Mv)) stop("corstr is assumed to be 'non_stat_M_dep' but Mv is not specified!")

if((corstr!="stat_M_dep" && corstr!="non_stat_M_dep") && !is.null(Mv))  stop("Mv is specified while corstr is assumed to be neither 
'stat_M_dep' nor 'non_stat_M_dep'!")

if(corstr=="non_stat_M_dep" && length(unique(nt)) !=1) stop("corstr cannot be assumed to be 'non_stat_M_dep' for unbalanced data!")

if(corstr=="unstructured" && length(unique(nt)) !=1) stop("corstr cannot be assumed to be 'unstructured' for unbalanced data!")

if(missing(R)) R<-NULL

if(corstr=="fixed" && is.null(R))  stop("corstr is assumed to be 'fixed' but R is not specified!")
if(corstr!="fixed" && !is.null(R)) stop("R is specified although corstr is not assumed to be 'fixed'!")

if(!is.null(R)) {
Rr <- nrow(R)
if(Rr != ncol(R)) stop("R is not square!")
if(Rr < maxclsz)  {stop("R is not big enough to accommodate some clusters!")} else
if(Rr > maxclsz)  {stop("R is larger than the maximum cluster!")}
}

if(missing(scale.fix))  scale.fix=FALSE
scale.fix <- as.integer(scale.fix)

if(missing(scale.value)) scale.value=1
scale.value<-as.integer(scale.value)

if(missing(eps)) eps=10^-6
eps<-as.double(eps)

if(missing(maxiter)) maxiter<-30
maxiter<-as.integer(maxiter)

if(missing(tol))  tol=10^-3
tol=as.double(tol)

if(missing(silent))  silent<-FALSE
silent<-as.integer(silent)

if (is.character(family)) family <- get(family)
if (is.function(family))  family <- family()

links <- c("identity","log","logit","inverse","probit","cloglog")
fams <- c("gaussian","poisson","binomial","Gamma","quasi")
varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")
corstrs <- c("independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", 
"AR-1", "unstructured")

linkv <- as.integer(match(c(family$link), links, -1))
if(linkv < 1) stop("unknown link!")

famv <- match(family$family, fams, -1)
if(famv < 1) stop("unknown family")
if(famv <= 4) varfunv <- famv
else varfunv <- match(family$varfun, varfuns, -1)
if(varfunv < 1) stop("unknown varfun!")

corstrv <- as.integer(match(corstr, corstrs, -1))
if(corstrv < 1) stop("unknown corstr!")

Mv <- as.integer(Mv)

if (!is.null(beta_int))
    {
        beta <- matrix(beta_int, nrow = 1)
        if(ncol(beta) != (K+1)) {stop("Dimension of beta != ncol(X)!")}
        message("user\'s initial regression estimate")
        
    }
    else {
        message("running glm to get initial regression estimate!")
### <tsl>	beta <- as.numeric(glm(m, family = family)$coef)
        mm <- match.call(expand.dots = FALSE)
        mm$R <- mm$beta_int <- mm$tol <- mm$maxiter <- mm$link <- 
        mm$varfun <-mm$corstr <- mm$Mv <- mm$silent <-mm$scale.fix <- 
        mm$scale.value <- mm$id<-
        mm$lambda <-mm$pindex<-mm$eps<-NULL
        mm[[1]]<-as.name("glm")
        beta <- eval(mm, parent.frame())$coef
### </tsl>
        print(beta)
       
}

beta_int=matrix(beta, ncol = 1)

beta_new<-beta_int

R.fi.hat=mycor_gee2(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R=R,scale.fix=scale.fix,scale.value=scale.fix)
Rhat=R.fi.hat$Ehat
fihat=R.fi.hat$fi

S.H.E.val=S_H_E_M(N,nt,y,X,K,family,beta_new,Rhat,fihat,lambda,pindex,eps)
S<-S.H.E.val$S
H<-S.H.E.val$H
E<-S.H.E.val$E

diff<-1
iter<-0

while(iter < maxiter) {

beta_old<-beta_new

beta_new<-matrix(beta_old)+ginv(H+N*E)%*%(S-N*E%*%matrix(beta_old))

R.fi.hat=mycor_gee2(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R,scale.fix,scale.value)
Rhat=R.fi.hat$Ehat
fihat=R.fi.hat$fi

S.H.E.M.val=S_H_E_M(N,nt,y,X,K,family,beta_new,Rhat,fihat,lambda,pindex,eps)
S<-S.H.E.M.val$S
H<-S.H.E.M.val$H
E<-S.H.E.M.val$E
M<-S.H.E.M.val$M

diff<-sum(abs(beta_old-beta_new)) 

iter<-iter+1
if (silent==1) cat("iter",iter,"beta_new",beta_new,"diff",diff,"\n")
if (diff <= tol) break
} #end of while

estb=beta_new
nv=naive.var<-ginv(H+N*E)
rv=robust.var<-ginv(H+N*E)%*%M%*%ginv(H+N*E)
final_iter=iter
final_diff=diff

fit <- list()
attr(fit, "class") <- c("PGee","gee","glm")
fit$title <- "PGEE: PENALIZED GENERALIZED ESTIMATING EQUATIONS FOR LONGITUDINAL DATA"
fit$version <- "Version: 1.1"
links <- c("Identity", "Logarithm", "Logit", "Reciprocal", "Probit","Cloglog")
varfuns <- c("Gaussian", "Poisson", "Binomial", "Gamma")
corstrs <- c("Independent", "Fixed", "Stationary M-dependent",
              "Non-Stationary M-dependent", "Exchangeable", "AR-1",
              "Unstructured")
fit$model <- list()
fit$model$link <- links[linkv]
fit$model$varfun <- varfuns[varfunv]
fit$model$corstr <- corstrs[corstrv]
if(!is.na(match(c(corstrv), c(3, 4))))
fit$model$M <- Mv
fit$call <- call
fit$terms <- Terms
fit$formula <- as.vector(attr(Terms, "formula"))
#fit$contrasts <- attr(X, "contrasts")
fit$nobs <- nobs
fit$iterations <- final_iter
fit$coefficients <- as.vector(estb)
fit$nas <- is.na(fit$coefficients)
names(fit$coefficients) <- xnames
eta <- as.vector(X %*% fit$coefficients)
fit$linear.predictors <- eta
##Rchange
mu <- as.vector(family$linkinv(eta))
##
fit$fitted.values <- mu
fit$residuals <- y - mu
fit$family <- family
fit$y <- as.vector(y)
fit$id <- as.vector(id)
fit$max.id <- maxcl
fit$working.correlation <- Rhat[1:maxclsz,1:maxclsz,which(avec==maxclsz)[1]]
fit$scale <- fihat
fit$epsilon<-eps
fit$lambda.value<-lambda
fit$robust.variance <- rv
fit$naive.variance <- nv
fit$xnames <- xnames
fit$error <- final_diff
dimnames(fit$robust.variance) <- list(xnames, xnames)
dimnames(fit$naive.variance) <- list(xnames, xnames)
fit

}
