`RM` <-
function(X, W, se = TRUE, sum0 = TRUE, etaStart)
{
#...X: 0/1 person*item matrix

#-------------------main programm-------------------

call<-match.call()
groupvec <- 1
mpoints <- 1
model <- "RM"

if (missing(W)) W <- NA
else W <- as.matrix(W)

if (missing(etaStart)) etaStart <- NA
else etaStart <- as.vector(etaStart)

XWcheck <- datcheck(X,W,mpoints,groupvec,model)                              #inital check of X and W
X <- XWcheck$X

lres <- likLR(X,W,mpoints,groupvec,model,st.err=se,sum0,etaStart)
parest <- lres$parest                             #full groups for parameter estimation

loglik <- -parest$minimum                         #log-likelihood value
iter <- parest$iterations                         #number of iterations
convergence <- parest$code
etapar <- parest$estimate                         #eta estimates
betapar <- as.vector(lres$W%*% etapar)            #beta estimates
if (se) {
  se.eta <- sqrt(diag(solve(parest$hessian)))         #standard errors
  se.beta <- sqrt(diag(lres$W%*%solve(parest$hessian)%*%t(lres$W)))   #se beta
} else {
  se.eta <- rep(NA,length(etapar))
  se.beta <- rep(NA,length(betapar))
}

X01 <- lres$X01
labs <- labeling.internal(model,X,X01,lres$W,etapar,betapar,mpoints,max(groupvec))    #labeling for L-models
W <- labs$W
etapar <- labs$etapar
betapar <- labs$betapar

etapar <- -etapar          # output difficulty  rh 25-03-2010

npar <- dim(lres$W)[2]                            #number of parameters

result <- list(X=X,X01=X01,model=model,loglik=loglik,npar=npar,iter=iter,convergence=convergence,
               etapar=etapar,se.eta=se.eta,hessian=parest$hessian,betapar=betapar,
               se.beta=se.beta,W=W,call=call)
class(result) <- c("dRm","Rm","eRm")                    #classes: dichotomous RM, RM (RM, PCM, RSM), and extended RM (all)
result
}

