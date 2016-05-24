
pause <- function() {
  cat("Press ENTER/RETURN/NEWLINE to continue.")
  readLines(n=1)
  invisible()
}

######################################################
#
# This is a fixed design simulation example for the GS/UPS.
# 
######################################################

### set up the parameters of the simulation example
p <- 1000
n <- p
### The parameters we used for simulation
v <- 0.4
r <- 3
### The true minimal signal strength
tau <- sqrt(2*r*log(p))
### The true sparsity level
sp <- p^(1-v)

pause()

### use the square root of Omega as the predictor matrix
### where Omega is a tri-diagonal matrix with off-diagonal element rho
rho <- 0.45
ii <- c(1:p,1:(p-1),2:p)
jj <- c(1:p,2:p,1:(p-1))
xx<-c(rep(1,p),rho*rep(1,2*(p-1)))
Omega <- sparseMatrix(ii,jj,x = xx) 
eigenOmega <- eigen(Omega)
OmegaVec <- eigenOmega$vectors
OmegaVal <- eigenOmega$values
OmegaRoot <- OmegaVec %*% diag(sqrt(OmegaVal)) %*% t(OmegaVec)

X <- OmegaRoot

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
delta <- 0.01/log(p)
gram.gram<-ThresholdGram(gram,delta)
gram<-gram.gram$gram
gram.bias<-gram.gram$gram.bias

pause()

### search for all the connected subgraphs with no more than nm nodes
nm <- 3
neighbor <- (gram!=0)
cg.all <- FindAllCG(neighbor,nm)

pause()

### GS for fixed design
y.tilde<-t(X)%*%Y

### The essential tuning parameters of the screening step (v, r) 
### are tied to the sparsity level and the minimal signal strength
### sp = p^(1-v) and tau = sqrt(2*r*log(p)).
### Here we use the (v, r) in the generative model from which 
### we generate the data. This is the oracle case.
survivor <- ScreeningStep(y.tilde, gram, cg.all, nm, v, r)

### lambda and uu are the tuning parameters of the cleaning step, 
### which are simple functions of the essential tuning parameters (v, r).
lambda <- sqrt(2*v*log(p))
uu <- tau
estimate.gs <- CleaningStep(survivor, y.tilde, gram, lambda, uu)

ham.gs <- sum(signbeta != sign(estimate.gs))

### the estimated sparse level    
sum(estimate.gs!=0)

### hamming error of the graphlet screening
ham.gs

pause()

### when nm = 1 in the screening step, it is identical to that of UPS.
survivor <- ScreeningStep(y.tilde, gram, cg.all, 1, v, r)
estimate.ups <- CleaningStep(survivor, y.tilde, gram, lambda, uu)
ham.ups <- sum(signbeta != sign(estimate.ups))

### the estimated sparse level of UPS   
sum(estimate.ups!=0)

### hamming error of UPS
ham.ups


pause()

### GS for fixed design, perform GS with perturbed parameters. 
### Now, instead of using the true parameters (v, r) in the generative model
### from which we simulate the data, we use perturbated (v, r) to tune the 
### graphlet screening. This example shows the robustness of graphlet screening
### against the moderate estimation error in tuning paramter estimation.

### In this example, the perturbed parameters (vp, rp) are obtained 
### by adding/substracting 10% to/from the real parameters (v, r).
vp <- v*(1+0.1*sample(c(-1,1),1))
rp <- r*(1+0.1*sample(c(-1,1),1))
vp
rp
survivor.perturb <- ScreeningStep(y.tilde, gram, cg.all, nm, vp, rp)

### The tuning parameters of the cleaning step are tied to (vp, rp)
lambdap <- sqrt(2*vp*log(p))
uup <- sqrt(2*rp*log(p))
estimate.gs.perturb <- CleaningStep(survivor.perturb, y.tilde, gram, lambdap, uup)

ham.gs.perturb <- sum(signbeta != sign(estimate.gs.perturb))

### the estimated sparse level    
sum(estimate.gs.perturb!=0)

### hamming error of the graphlet screening with perturbed parameters.
ham.gs.perturb
