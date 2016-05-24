cat("\ntest CAR and SEM:")

data(scotlip)

corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
          data=scotlip,family=poisson(),
          adjMatrix=Nmatrix) ## 4 s.
## same without optim: run in scotlip examples; cf also autoregressive.Rd for ML fits

set.seed(124)
ldl <- selfAdjointSolverCpp(Nmatrix)
Lmat <- ldl$u %*% diag(sqrt(1/(1-0.17*ldl$d)))
lp <- 0.1 + 3* Lmat %*% rnorm(ncol(Lmat)) ## single intercept beta =0.1; lambda=3
resp <- rbinom(ncol(Lmat),1,1/(1+exp(-lp)))
donn <- data.frame(npos=resp,nneg=1-resp,gridcode=scotlip$gridcode)

# CAR by Laplace with 'outer' estimation of rho
corrHLfit(cbind(npos,nneg)~1 +adjacency(1|gridcode),
          adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="ML") ## 43 s.

# CAR by Laplace with 'inner' estimation of rho
HLCor(cbind(npos,nneg)~1 +adjacency(1|gridcode),
          adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="ML") ## 4 s.

# CAR by SEMs +optimsmooth ... slow
#corrHLfit(cbind(npos,nneg)~1 +adjacency(1|gridcode),
#          adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="SEM")

# CAR by single SEM
HLCor(cbind(npos,nneg)~1 +adjacency(1|gridcode),
          adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="SEM") ## 5 s.
