gammahetero1 <-  function(formula1,formula2){

# Mean model matrix construction
tr=terms(formula1)
data.1 = eval(attr(tr,"variables"))
Y <- data.1[[1]]
X <- model.matrix(formula1)
k <- ncol(X)
n <- nrow(X)

# Dispersion model matrix construction
Z <- model.matrix(formula2)
r <- ncol(Z)


# Functions for mu, alpha, W_g and W_b

mu <- function(betas){ exp(X%*%betas)}
alpha <- function(gammas){ exp(Z%*%gammas) }
Sigma.b <- function(betas,gammas){
 alpha(gammas)
}
Sigma.g <- function(gammas){
 (alpha(gammas)^(-2))*(trigamma(alpha(gammas))-alpha(gammas)^(-1))^(-1)
}

Y.tilde <- function(betas){
Ytilde = X%*%betas + Y/mu(betas) - 1
return(Ytilde)
}


Y.mono <- function(betas,gammas){
mu = mu(betas)
alpha = alpha(gammas)
Ymono = Z%*%gammas - (1/alpha)*((trigamma(alpha)-(1/alpha))^(-1))*(digamma(alpha)-
log(alpha*Y/mu)- 1 + Y/mu)
return(Ymono)
}


# Initial Values for Beta's
Y.trans = exp(Y)
beta <- solve(t(X)%*%X) %*% t(X)%*%Y.trans

# Initial Values for Gamma's
Y.hat <- X%*%beta
error <- (Y.hat - Y)^2
glm <- glm.fit(Z,error, family=Gamma(log))
gamma <- glm$coefficients


# Iterative Procedure
 convergencia <- 1
 iteracion <- 0
 while(convergencia>0.00001){
  beta.m = beta
  gamma.m = gamma
  Sigmab = diag(as.vector(Sigma.b(beta.m,gamma.m)))
  Ytilde = Y.tilde(beta.m)
  beta = (solve(t(X)%*%Sigmab%*%X))%*%t(X)%*%Sigmab%*%Ytilde
  Sigmag = diag(as.vector(1/Sigma.g(gamma.m)))
  Ymono = Y.mono(beta,gamma.m)
  gamma = solve(t(Z)%*%Sigmag%*%Z) %*% (t(Z)%*%Sigmag%*%Ymono)
  convergencia = sum(((gamma-gamma.m)/gamma.m)^2)
  iteracion = 1+iteracion
 }

# Variance and Covariance Matrix for Mean Model
W.b = diag(as.vector(alpha(gamma)))
cov.betas = solve(t(X)%*%W.b%*%X)
# Information Matrix for Dispersion Model
W.g = diag(as.vector((1-trigamma(alpha(gamma))-1/alpha(gamma)^2)*exp(2*Z%*%gamma)))
cov.gammas = solve( t(Z)%*%W.g%*%Z )

# AIC

lik <- function(mu,alpha){
 sum(alpha*log(alpha/mu) - lgamma(alpha) - (alpha-1)*log(Y)-Y*mu/alpha)
}

AIC = -2*lik(mu(beta),alpha(gamma)) + 2*(k+r) 

ICL = beta - qt(1-0.05/2,df=n-r-1)*sqrt(diag(cov.betas))
ICR = beta + qt(1-0.05/2,df=n-r-1)*sqrt(diag(cov.betas))
ICB = cbind(ICL,ICR)

GCL = gamma - qt(1-0.05/2,df=n-r-1)*sqrt(diag(cov.gammas))
GCR = gamma + qt(1-0.05/2,df=n-r-1)*sqrt(diag(cov.gammas))
ICG=cbind(GCL,GCR)
 
#Outputs
list(X=X,Z=Z,beta=beta,gamma=gamma,ICB=ICB, ICG=ICG,
CovarianceMatrixbeta=cov.betas,CovarianceMatrixgamma =cov.gammas,
AIC=AIC,iteration=iteracion,convergence=convergencia)
}


