
pause <- function() {
  cat("Press ENTER/RETURN/NEWLINE to continue.")
  readLines(n=1)
  invisible()
}

######################################################
#
# This is a random design simulation example for the GS/UPS
# 
######################################################

### set up the parameters of the simulation example
p <- 3000
theta <- 0.975
n <- round(p^theta)
v <- 0.35
r <- 3.5
tau <- sqrt(2*r*log(p))
sp <- p^(1-v)

pause()

### generate the predictor matrix from N(0,Omega/n)
### where Omega is a tri-diagonal matrix with off-diagonal element rho
rho <- 0.4
ii <- c(1:p,1:(p-1),2:p)
jj <- c(1:p,2:p,1:(p-1))
xx<-c(rep(1,p),rho*rep(1,2*(p-1)))
Omega <- sparseMatrix(ii,jj,x = xx) 
eigenOmega <- eigen(Omega)
OmegaVec <- eigenOmega$vectors
OmegaVal <- eigenOmega$values
OmegaRoot <- OmegaVec %*% diag(sqrt(OmegaVal)) %*% t(OmegaVec)

X <- matrix(rnorm(n*p),n,p)
X <- X %*% OmegaRoot/sqrt(n)

pause()

### generate the signals as in Jin and et al (2012) 
uun <- rnorm(p)^2 * (runif(p)<0.2)/6
taubeta <- (1+uun)*tau # /diag.scale
signbeta <- sample(c(-1,1),p,replace=T)
supportbeta <- (runif(p)<p^(-v))
signbeta<-supportbeta*signbeta
beta <- signbeta*taubeta

### number of true signals
sum(supportbeta)

pause()

### generate the response variable
noise <- rnorm(n)
Y <- X %*% beta + noise

pause()

### threshold the gram matrix 
gram<-t(X)%*%X
delta <- 1/log(p)
gram.gram<-ThresholdGram(gram,delta)
gram<-gram.gram$gram
gram.bias<-gram.gram$gram.bias

pause()

### search for all the connected subgraphs with no more than nm nodes
nm <- 3
neighbor <- (gram!=0)
cg.all <- FindAllCG(neighbor,nm)

pause()

###iterative GS for random design
y.tilde<-t(X)%*%Y
### iterative GS use the sparse level sp and the minimal signals strength tau
### to calculate the tuning parameters of the screening step and the cleaning step.
### More specifically, v = 1-log(sp)/log(p), and r = tau^2/2/log(p)

### In this example, we use the true sparse level and the true minimal signal strength
### of the generative model from which the data is generated. This is the oracle case.
result.IterGS<-IterGS(y.tilde,gram,gram.bias,cg.all, sp,tau, nm)
estimate.gs<-result.IterGS$estimate
result.IterGS$n.iter

ham.gs <- sum(signbeta != sign(estimate.gs))

### the estimated sparse level    
sum(estimate.gs!=0)

### hamming error of the graphlet screening
ham.gs


pause()

### iterative GS for random design, repeat the above example with perturbed parameters. 
###  sparse level and the minimal signal strength, we use perturbated (sp, tau) to tune
### graphlet screening. This example shows the robustness of graphlet screening
### against the moderate estimation error in tuning paramter estimation.

### In this example, the perturbed parameters (sp.perturb, tau.perturb) are obtained 
### by adding/substracting 10% to/from the real parameters (sp, tau). 
sp.perturb <- sp*(1+0.1*sample(c(-1,1),1))
tau.perturb <- tau*(1+0.1*sample(c(-1,1),1))
sp.perturb
tau.perturb

result.IterGS.perturb<-IterGS(y.tilde,gram,gram.bias,cg.all, sp.perturb,tau.perturb, nm)
estimate.gs.perturb<-result.IterGS.perturb$estimate
result.IterGS.perturb$n.iter
ham.gs.perturb <- sum(signbeta != sign(estimate.gs.perturb))

### the estimated sparse level    
sum(estimate.gs.perturb!=0)

### hamming error of the graphlet screening with perturbed tuning parameters
ham.gs.perturb


