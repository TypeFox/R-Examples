JGee1 <-
function(formula,id,data,nr,na.action=NULL,family=gaussian(link="identity"),
corstr1="independence",Mv=NULL,corstr2="independence",
beta_int=NULL,R1=NULL,R2=NULL,scale.fix=FALSE,
scale.value=1,maxiter=25,tol=10^-3,silent=FALSE)  {

call <- match.call()
m <- match.call(expand.dots = FALSE)

m$family <-m$link <- m$varfun <-
m$nr<-m$beta_int <-
m$corstr1 <- m$Mv<- m$corstr2 <-m$R1 <- m$R2 <-
m$scale.fix <- m$scale.value <-
m$maxiter <- m$tol <-m$silent <-NULL

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

id<-model.extract(m, id)

if(is.null(id)) {
stop("Id variable not found!")
}

if(is.null(nr)) {
stop("nr variable not found!")
}

if(length(id) != length(y))  stop("Id and y not same length!")

if(!(is.double(X)))  X <- as.double(X)
if(!(is.double(y)))  y <- as.double(y)
if(!(is.double(id))) id <- as.double(id)

N<-length(unique(id))
nr<-as.integer(nr)
K<-ncol(X)-1

xnames <- dimnames(X)[[2]]
if(is.null(xnames)) {
xnames <- paste("x", 0:K, sep = "")
dimnames(X) <- list(NULL, xnames)
}


avec <- as.integer(unlist(lapply(split(id, id), "length")))/nr
maxclsz <-max(avec)
maxcl <- maxclsz
nt<-avec
nobs<-sum(nt)*nr

if(missing(family)) family=gaussian(link="identity")

if(missing(corstr1)) corstr1="independence"

if(missing(corstr2)) corstr2="independence"

if(missing(Mv)) Mv<-NULL

if(corstr1=="stat_M_dep" && is.null(Mv)) stop("corstr1 is assumed to be 'stat_M_dep' but Mv is not specified!")

if(corstr1=="non_stat_M_dep" && is.null(Mv)) stop("corstr1 is assumed to be 'non_stat_M_dep' but Mv is not specified!")

if((corstr1!="stat_M_dep" && corstr1!="non_stat_M_dep") && !is.null(Mv)) stop("Mv is specified while corstr1 is assumed to be neither 
'stat_M_dep' nor 'non_stat_M_dep'!")

if(corstr1=="non_stat_M_dep" && length(unique(nt)) !=1) stop("corstr1 cannot be assumed to be 'non_stat_M_dep' for unbalanced data!")

if(corstr1=="unstructured" && length(unique(nt)) !=1) stop("corstr1 cannot be assumed to be 'unstructured' for unbalanced data!")

if(missing(R1)) R1<-NULL
if(missing(R2)) R2<-NULL

if(corstr1!="fixed" && corstr2=="fixed")  stop("corstr1 should also be assumed to be 'fixed' when corstr2='fixed'!")
if(corstr1=="fixed" && corstr2!="fixed")  stop("corstr2 should also be assumed to be 'fixed' when corstr1='fixed'!")

if(corstr1=="fixed" && corstr2=="fixed" && is.null(R1) && is.null(R2))   stop("corstr1 and corstr2 are assumed to be 'fixed' but R1 and R2 are not specified!")
if(corstr1=="fixed" && corstr2=="fixed" && is.null(R1) && !is.null(R2))  stop("corstr1 and corstr2 are assumed to be 'fixed' but R1 is not specified!")
if(corstr1=="fixed" && corstr2=="fixed" && !is.null(R1) && is.null(R2))  stop("corstr1 and corstr2 are assumed to be 'fixed' but R2 is not specified!")

if(corstr1!="fixed" && corstr2!="fixed" && !is.null(R1)&& !is.null(R2)) stop("R1 and R2 are specified although corstr1 and corstr2 are not assumed to be 'fixed'!")
if(corstr1!="fixed" && corstr2!="fixed" && !is.null(R1)&& is.null(R2))  stop("R1 is specified although corstr1 and corstr2 are not assumed to be 'fixed'!")
if(corstr1!="fixed" && corstr2!="fixed" && is.null(R1)&& !is.null(R2))  stop("R2 is specified although corstr1 and corstr2 are not assumed to be 'fixed'!")

if(!is.null(R1)) {
Rr1 <- nrow(R1)
if(Rr1 != ncol(R1)) {stop("R1 is not square!")} else
if(Rr1 < maxclsz)   {stop("R1 is not big enough to accommodate some clusters!")} else
if(Rr1 > maxclsz)   {stop("R1 is larger than the maximum cluster!")}
}

if(!is.null(R2)) {
Rr2 <- nrow(R2)
if(Rr2 != ncol(R2)) {stop("R2 is not square!")} else
if(Rr2 < nr)  {stop("R2 is not big enough to accommodate some clusters!")} else
if(Rr2 > nr)  {stop("R2 is larger than the number of responses!")} 
}

if(missing(scale.fix))  scale.fix=FALSE
scale.fix <- as.integer(scale.fix)

if(missing(scale.value)) scale.value=1
scale.value<-as.integer(scale.value)

if(missing(maxiter)) maxiter<-25
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
corstrs1 <- c("independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", 
"AR-1", "unstructured")
corstrs2 <- c("independence", "fixed", "exchangeable","unstructured")

linkv <- as.integer(match(c(family$link), links, -1))
if(linkv < 1) stop("unknown link!")

famv <- match(family$family, fams, -1)
if(famv < 1) stop("unknown family")
if(famv <= 4) varfunv <- famv
else varfunv <- match(family$varfun, varfuns, -1)
if(varfunv < 1) stop("unknown varfun!")

corstrv1 <- as.integer(match(corstr1, corstrs1, -1))
if(corstrv1 < 1) stop("unknown corstr1!")

corstrv2 <- as.integer(match(corstr2, corstrs2, -1))
if(corstrv2 < 1) stop("unknown corstr2!")

Mv <- as.integer(Mv)

if (!is.null(beta_int))
    {
        beta <- matrix(beta_int, nrow = 1)
        if(ncol(beta) != (K+1)) {stop("Dimension of beta != ncol(X)!")}
        message("user\'s initial regression estimate")
        
    }
    else {
        message("running gee to get initial regression estimate")
### <tsl>	beta <- as.numeric(gee(m,family = family,corstr="independence")$coefficients)
        mm <- match.call(expand.dots = FALSE)
        mm$nr<- mm$R1 <-mm$R2 <- mm$beta_int <- mm$tol <- mm$maxiter <- mm$link <- mm$varfun <-mm$corstr1 <- 
        mm$corstr2 <-mm$Mv <- mm$silent <-mm$scale.fix <- mm$scale.value <- mm$id<- NULL
        mm[[1]]<-as.name("glm")
        beta <- eval(mm, parent.frame())$coefficients
### </tsl>
        print(beta)
       
}

beta_int=matrix(beta, ncol = 1)

beta_new<-beta_int

R.fi.hat=mycor_jgee1(N,nr,nt,y,X,family,beta_new,corstr1,Mv,corstr2,maxclsz,R1=R1,R2=R2,scale.fix=scale.fix,scale.value=scale.fix)
Rhat=R.fi.hat$Ehat
fihat=R.fi.hat$fi

S.H.E.val=S_H1(N,nr,nt,y,X,K,family,beta_new,Rhat,fihat)
S<-S.H.E.val$S
H<-S.H.E.val$H

diff<-1
iter<-0

while(iter < maxiter) {

beta_old<-beta_new

beta_new<-matrix(beta_old)+(ginv(H)%*%S)

R.fi.hat=mycor_jgee1(N,nr,nt,y,X,family,beta_new,corstr1,Mv,corstr2,maxclsz,R1,R2,scale.fix,scale.value)
Rhat=R.fi.hat$Ehat
fihat=R.fi.hat$fi
C1hat=R.fi.hat$cor1
C2hat=R.fi.hat$cor2

S.H.E.M.val=S_H1(N,nr,nt,y,X,K,family,beta_new,Rhat,fihat)
S<-S.H.E.M.val$S
H<-S.H.E.M.val$H
M<-S.H.E.M.val$M

diff<-sum(abs(beta_old-beta_new)) 

iter<-iter+1
if (silent==1) cat("iter",iter,"diff",diff,"\n")
if (diff <= tol) break
} #end of while

estb=beta_new
nv=naive.var<-ginv(H)
rv=robust.var<-ginv(H)%*%M%*%ginv(H)
final_iter=iter
final_diff=diff

fit <- list()
attr(fit, "class") <- c("JGee1","gee","glm")
fit$title <- "JGEE: JOINT GENERALIZED ESTIMATING EQUATIONS FOR CLUSTERED DATA"
fit$version <- "Version: 1.0"
links <- c("Identity", "Logarithm", "Logit", "Reciprocal", "Probit","Cloglog")
varfuns <- c("Gaussian", "Poisson", "Binomial", "Gamma")
corstrs1 <- c("Independent", "Fixed", "Stationary M-dependent",
              "Non-Stationary M-dependent", "Exchangeable", "AR-1",
              "Unstructured")
corstrs2 <- c("Independent", "Fixed","Exchangeable", "Unstructured")
fit$model <- list()
fit$model$link <- links[linkv]
fit$model$varfun <- varfuns[varfunv]
fit$model$corstr1 <- corstrs1[corstrv1]
if(!is.na(match(c(corstrv1), c(3, 4))))
fit$model$M <- Mv
fit$model$corstr2 <- corstrs2[corstrv2]
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
fit$working.correlation1 <- C1hat[1:maxclsz,1:maxclsz,which(avec==maxclsz)[1]]
fit$working.correlation2 <- C2hat
fit$working.correlation <- Rhat[1:(maxclsz*nr),1:(maxclsz*nr),which(avec==maxclsz)[1]]
fit$scale <- fihat
fit$robust.variance <- rv
fit$naive.variance <- nv
fit$xnames <- xnames
fit$error <- final_diff
dimnames(fit$robust.variance) <- list(xnames, xnames)
dimnames(fit$naive.variance) <- list(xnames, xnames)
fit

}
