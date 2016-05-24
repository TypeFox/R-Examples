mcmc.fmem <-
function(params){

rmvnorm.l <- function(mean, sigma){
d <- length(mean)
sigma2 <- chol(sigma)
t(mean + t(sigma2)%*%rnorm(d))
}


Wishart <- function(df,Omega){
  p <- ncol(Omega)
  x <- matrix(rnorm(df*p),df,p)
  chole <- chol(Omega)
  t(chole)%*%(t(x)%*%x)%*%chole
}  

burn.in <- params$burn.in
post.sam.s <- params$post.sam.s
thin <- params$thin
y <- params$y
p <- params$p
q <- params$q
ks <- params$ks
u <- params$u
pdf <- params$pdf
cdf <- params$cdf
n <- params$n
kappa <- params$kappa
omeg <- 1/params$omeg
family <- params$family
homo <- params$homo
if(homo == 0) heter <- params$heter
 
total <- burn.in + post.sam.s*thin
ancho <- floor(seq(2, total, length=10))

bar <- txtProgressBar(min=0, max=ancho[10], initial=0, width=50, char="+", style=3)

mu_m0 <- rep(0,q)
sigma2_mu0 <- diag(q)*1000
sigma2_mu0.I <- solve(sigma2_mu0)

nu_m <- q
omega_m <- diag(q)*1000

a.sigma2_y <- 0.0001/2
b.sigma2_y <- 0.0001/2

rho.a <- matrix(0,total,q)
rho0 <- rep(0,q)
S.rho <- diag(q)*1000
S.rho.I <- solve(S.rho)
rho.a[1,] <- params$rho.i
M <- params$M
params$rho0 <- rho0
params$S.rho <- S.rho

sigma2_m.a <- array(0,c(q,q,total))
sigma2_m.a2 <- array(0,c(q,q,total))
mu_m.a <- matrix(0,total,q)
m_i.a <- matrix(0,n,q)
sigma2_y.a <- matrix(0,total,1)

sigma2_m.a[,,1] <- as.matrix(var(M))
mu_m.a[1,] <- apply(M,2,mean)
m_i.a <- M
sigma2_y.a[1] <- params$sigma2_y

beta.a <- matrix(0,total,p)
beta0 <- rep(0,p)
S.beta <- diag(p)*1000
S.beta.II <- matrix(0,p+q,p+q)
S.beta.II[1:p,1:p] <- solve(S.beta)
S.beta.II[(p+1):(p+q),(p+1):(p+q)] <- S.rho.I
beta_rho0 <- matrix(0,p+q,1)
beta_rho0[1:p] <- beta0
beta_rho0[(p+1):(q+p)] <- rho0

beta.a[1,] <- params$beta.i

X <- params$X
params$beta0 <- beta0
params$S.beta <- S.beta

tau.a <- matrix(0,total,length(ks))
a.tau <- 0.001
b.tau <- 0.001
if(sum(ks) > 0){
   alpha.a <- matrix(0,total,sum(ks))
   alpha.a[1,] <- params$alpha.i
   B <- params$B
   }
else{alpha.a <- matrix(0,total,1)
       B <- matrix(1,n,1)
   }

cont <- 1
uii <- matrix(0,n,1)
inv.omega_m <- solve(omega_m)


if(family!="Normal" && family!="Laplace" && attr(params$eta,"know")==0){
extra.parameter <- params$extra.parameter
nu0 <- params$nu0
nu.a <- matrix(0,total,length(nu0))
    nu.a[1,] <- nu0
}else{
nu.a <- kronecker(matrix(1,total,length(params$eta)),t(params$eta))
}


for(l in 2:total){

M_m <- 0
omega_m.a <- matrix(0,ncol(M),ncol(M))

if(homo == 1){
ss <- (y - X%*%beta.a[l-1,] - m_i.a%*%rho.a[l-1,]- B%*%alpha.a[l-1,])^2/sigma2_y.a[l-1] + apply((M - m_i.a)^2,1,sum)/(sigma2_y.a[l-1]*omeg) +
      diag((m_i.a-kronecker(t(mu_m.a[l-1,]),matrix(1,n,1)))%*%tcrossprod(sigma2_m.a[,,l-1],(m_i.a-kronecker(t(mu_m.a[l-1,]),matrix(1,n,1)))))
}
else{
ss <- (y - X%*%beta.a[l-1,] - m_i.a%*%rho.a[l-1,]- B%*%alpha.a[l-1,])^2/heter$sigma2y + apply((M - m_i.a)^2/(heter$sigma2xi),1,sum) +
      diag((m_i.a-kronecker(t(mu_m.a[l-1,]),matrix(1,n,1)))%*%tcrossprod(sigma2_m.a[,,l-1],(m_i.a-kronecker(t(mu_m.a[l-1,]),matrix(1,n,1)))))
}

if(family!="Hyperbolic" && family!="Laplace") uii <- u(ss,nu.a[l-1,])

if(homo == 1){
sig.m_i.a1 <- solve((omeg*sigma2_y.a[l-1])^(-1)*diag(ncol(M)) + sigma2_m.a[,,l-1]  + rho.a[l-1,]%*%t(rho.a[l-1,])/(sigma2_y.a[l-1]))
chole <- chol(sig.m_i.a1)
}

for(i in 1:n){
if(family=="Hyperbolic" || family=="Laplace") uii[i] <- u(ss[i],nu.a[l-1,])

if(homo == 1){
sig.m_i.a <- kappa(uii[i])*sig.m_i.a1

mu.m_i.a <- sig.m_i.a%*%(sigma2_m.a[,,l-1]%*%mu_m.a[l-1,]/kappa(uii[i]) + (kappa(uii[i])*omeg*sigma2_y.a[l-1])^(-1)*diag(ncol(M))%*%M[i,]
+ rho.a[l-1,]*(y[i]-X[i,]%*%beta.a[l-1,] - B[i,]%*%alpha.a[l-1,])/(sigma2_y.a[l-1]*kappa(uii[i])))
m_i.a[i,] <- mu.m_i.a + sqrt(kappa(uii[i]))*crossprod(chole,rnorm(q))

}else{
dfg <- diag(q)
diag(dfg) <- 1/heter$sigma2xi[i,]
sig.m_i.a <- kappa(uii[i])*solve( dfg + sigma2_m.a[,,l-1]  + rho.a[l-1,]%*%t(rho.a[l-1,])/(heter$sigma2y[i]))
chole <- chol(sig.m_i.a)

mu.m_i.a <- sig.m_i.a%*%(sigma2_m.a[,,l-1]%*%mu_m.a[l-1,]/kappa(uii[i]) + M[i,]/(heter$sigma2xi[i,]*kappa(uii[i]))
+ rho.a[l-1,]*(y[i]-X[i,]%*%beta.a[l-1,] - B[i,]%*%alpha.a[l-1,])/(heter$sigma2y[i]*kappa(uii[i])))

m_i.a[i,] <- mu.m_i.a + crossprod(chole,rnorm(q))
}

omega_m.a <- omega_m.a + tcrossprod((m_i.a[i,]-mu_m.a[l-1,]),(m_i.a[i,]-mu_m.a[l-1,]))/kappa(uii[i])
M_m <- M_m + crossprod((M[i,]- m_i.a[i,]),(M[i,] - m_i.a[i,]))/kappa(uii[i])
}


if(family!="Normal" && family!="Laplace"  && attr(params$eta,"know")==0){
nu.a[l,] <- extra.parameter(nu.a[l-1,], uii, ss) 
}

if(homo ==1){
MX <- cbind(X,m_i.a)
sig.beta.a <- solve(S.beta.II + crossprod(MX,matrix(1/(sigma2_y.a[l-1] * kappa(uii)),n,(p+q))*MX))
mu.beta.a <-  sig.beta.a%*%(S.beta.II%*%beta_rho0 + crossprod(MX,matrix(1/(sigma2_y.a[l-1]*kappa(uii)),n,1)*(y-B%*%alpha.a[l-1,])))
temp <- rmvnorm.l(mean=mu.beta.a, sigma=sig.beta.a)
beta.a[l,] <- temp[1:p]
rho.a[l,] <- temp[(p+1):(q+p)]
if(sum(ks) > 0){
alpha.a[l,] <- alpha.a[l-1,]

if(length(ks)==1){
alpha0 <- rep(0,ks)
scale.tau <- 2*b.tau + crossprod(alpha.a[l-1,],alpha.a[l-1,])
tau.a[l,1] <- 1/rgamma(1, shape= (ks/2 + a.tau), scale= 2/scale.tau)

sig.alpha.a <- solve(diag(ks)/tau.a[l,1] + crossprod(B,matrix(1/(sigma2_y.a[l-1]*kappa(uii)),n,ks)*B))
mu.alpha.a <- sig.alpha.a%*%(crossprod(B,(1/(sigma2_y.a[l-1]*kappa(uii)))*(y-X%*%beta.a[l,]-m_i.a%*%rho.a[l,] )))
alpha.a[l,] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)
}
else{
lmi <- 1
for(i in 1:length(ks)){
alpha0 <- rep(0,ks[i])
    B.i <- B[,lmi:(lmi+ks[i]-1)] 
scale.tau <- 2*b.tau + crossprod(alpha.a[l-1,lmi:(lmi+ks[i]-1)],alpha.a[l-1,lmi:(lmi+ks[i]-1)])
tau.a[l,i] <- 1/rgamma(1, shape= (ks[i]/2 + a.tau), scale= 2/scale.tau)

sig.alpha.a <- solve(diag(ks[i])/tau.a[l,i] + crossprod(B.i,matrix(1/(sigma2_y.a[l-1]*kappa(uii)),n,ks[i])*B.i))
mu.alpha.a <- sig.alpha.a%*%(crossprod(B.i,(1/(sigma2_y.a[l-1]*kappa(uii)))*(y-X%*%beta.a[l,]-m_i.a%*%rho.a[l,] - B[,-c(lmi:(lmi+ks[i]-1))]%*%alpha.a[l,-c(lmi:(lmi+ks[i]-1))])))
alpha.a[l,lmi:(lmi+ks[i]-1)] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)

lmi <- lmi + ks[i]
}
}
}else{tau.a[l] <- 0
          alpha.a[l,] <- 0}


}else{
MX <- cbind(X,m_i.a)
sig.beta.a <- solve(S.beta.II + crossprod(MX,matrix(1/(heter$sigma2y*kappa(uii)),n,(p+q))*MX))
mu.beta.a <-  sig.beta.a%*%(S.beta.II%*%beta_rho0 + crossprod(MX,matrix(1/(heter$sigma2y*kappa(uii)),n,1)*(y-B%*%alpha.a[l-1,])))
temp <- rmvnorm.l(mean=mu.beta.a, sigma=sig.beta.a)
beta.a[l,] <- temp[1:p]
rho.a[l,] <- temp[(p+1):(q+p)]

if(sum(ks) > 0){
alpha.a[l,] <- alpha.a[l-1,]

if(length(ks)==1){
alpha0 <- rep(0,ks)
scale.tau <- 2*b.tau + crossprod(alpha.a[l-1,],alpha.a[l-1,])
tau.a[l,1] <- 1/rgamma(1, shape= (ks/2 + a.tau), scale= 2/scale.tau)

sig.alpha.a <- solve(diag(ks)/tau.a[l,1] + crossprod(B,matrix(1/(heter$sigma2y*kappa(uii)),n,ks)*B))
mu.alpha.a <- sig.alpha.a%*%(crossprod(B,(1/(heter$sigma2y*kappa(uii)))*(y-X%*%beta.a[l,]-m_i.a%*%rho.a[l,] )))
alpha.a[l,] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)
}
else{
lmi <- 1
for(i in 1:length(ks)){
alpha0 <- rep(0,ks[i])
    B.i <- B[,lmi:(lmi+ks[i]-1)] 
scale.tau <- 2*b.tau + crossprod(alpha.a[l-1,lmi:(lmi+ks[i]-1)],alpha.a[l-1,lmi:(lmi+ks[i]-1)])
tau.a[l,i] <- 1/rgamma(1, shape= (ks[i]/2 + a.tau), scale= 2/scale.tau)

sig.alpha.a <- solve(diag(ks[i])/tau.a[l,i] + crossprod(B.i,matrix(1/(heter$sigma2y*kappa(uii)),n,ks[i])*B.i))
mu.alpha.a <- sig.alpha.a%*%(crossprod(B.i,(1/(heter$sigma2y*kappa(uii)))*(y-X%*%beta.a[l,]-m_i.a%*%rho.a[l,] - B[,-c(lmi:(lmi+ks[i]-1))]%*%alpha.a[l,-c(lmi:(lmi+ks[i]-1))])))
alpha.a[l,lmi:(lmi+ks[i]-1)] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)

lmi <- lmi + ks[i]
}
}
}else{tau.a[l] <- 0
     alpha.a[l,] <- 0}
}

sig.mu_m.a <- solve(sigma2_mu0.I + sigma2_m.a[,,l-1]*n*mean(1/kappa(uii)))
mu.mu_m.a <-  sig.mu_m.a%*%(sigma2_mu0.I%*%mu_m0 + n*(sigma2_m.a[,,l-1])%*%apply(as.matrix(m_i.a)/matrix(kappa(uii),n,q),2,mean))
mu_m.a[l,] <- rmvnorm.l(mean=mu.mu_m.a, sigma=sig.mu_m.a)

sigma2_m.a[,,l] <- Wishart(n+q, solve(inv.omega_m + omega_m.a))
sigma2_m.a2[,,l] <- solve(sigma2_m.a[,,l])

if(homo==1){
scale.sigma2_y <- sum((y-X%*%beta.a[l,]-m_i.a%*%rho.a[l,]-B%*%alpha.a[l,])^2/kappa(uii)) + (1/omeg)*M_m + b.sigma2_y
sigma2_y.a[l] <- 1/rgamma(1,shape= ((n*(1+q) + a.sigma2_y)/2), scale=2/scale.sigma2_y)
}


if(l==ancho[cont]){
 Sys.sleep(0.5);
             setTxtProgressBar(bar,ancho[cont])
 cont <- cont + 1
}

l <- l + 1

}

close(bar)

size <- seq(burn.in+thin,total,length=post.sam.s)

if(attr(params$eta,"know")==0){
aa <- matrix(0,post.sam.s,p+3*q+sum(ks)+1+length(nu0)+length(ks))}
else{aa <- matrix(0,post.sam.s,p+3*q+sum(ks)+1+length(ks))}

cad <- as.vector(" ")

aa[,1:p] <- beta.a[size,]
cad <- cbind(cad,t(paste("beta",1:p)))

aa[,(p+1):(p+q)] <- rho.a[size,]
cad <- cbind(cad,t(paste("rho",1:q)))

aa[,(p+q+1):(p+2*q)] <- mu_m.a[size,]
cad <- cbind(cad,t(paste("mu_m",1:q)))

for(i in 1:q){
aa[,(p+2*q+i)] <- as.vector(sigma2_m.a2[i,i,size])
}
cad <- cbind(cad,t(paste("sigma2_m",1:q)))

if(homo==1){
aa[,(p+3*q+1)] <- sigma2_y.a[size]
cad <- cbind(cad,t(paste("sigma2_y")))
}

if(sum(ks) > 0){
aa[,(p+3*q+2):(p+3*q+1+sum(ks))] <- alpha.a[size,]
aa[,(p+3*q+sum(ks)+2):(p+3*q+sum(ks)+1+length(ks))] <- tau.a[size,]
cad <- cbind(cad,t(colnames(B)))
}


if(attr(params$eta,"know")==0){
if(ncol(nu.a)==1) {
cond <- var(nu.a[size])
if(cond==0) nu.a[size[1]] <- 1.1*mean(nu.a[size])
}
if(ncol(nu.a)==2){
cond <- diag(var(nu.a[size,]))
if(cond[1]==0) nu.a[size[1],1] <- 0.98*mean(nu.a[size,1])
if(cond[2]==0) nu.a[size[1],2] <- 0.98*mean(nu.a[size,2])
}
aa[,(p+3*q+sum(ks)+1+length(ks)+1):(p+3*q+sum(ks)+1+length(ks)+length(nu0))] <- nu.a[size,]
}


if(sum(ks)>0) cad <- cbind(cad,t(paste("tau_alpha",1:length(ks))))

if(family!="Normal" && family!="Laplace"  && attr(params$eta,"know")==0){
cad <- cbind(cad,t(paste("eta",1:length(nu0))))
}

haver <- apply(aa,2,var)
aa2 <- aa[,haver!=0]

cad <- cad[2:length(cad)]

colnames(aa2) <- cad

mi <- matrix(0,n,post.sam.s)

D_bar <- 0
sigma_mp <- matrix(0,q+1,q+1)

for(i in 1:post.sam.s){
sigma_mp[1,2:(q+1)] <- rho.a[(burn.in+i),]%*%sigma2_m.a2[,,(burn.in+i)]
sigma_mp[2:(q+1),1] <- t(sigma_mp[1,2:(q+1)])

if(homo==1){
sigma_mp[1,1] <- sigma2_y.a[(burn.in+i)] + crossprod(rho.a[(burn.in+i),],sigma2_m.a2[,,(burn.in+i)]%*%rho.a[(burn.in+i),])
sigma_mp[2:(q+1),2:(q+1)] <- omeg*sigma2_y.a[(burn.in+i)] + sigma2_m.a2[,,(burn.in+i)]
inv.sigma_mp <- solve(sigma_mp)
det_sigma_mp <- det(sigma_mp)
}

comp_s <- X%*%beta.a[(burn.in+i),] + B%*%alpha.a[(burn.in+i),] + sum(rho.a[(burn.in+i),]*mu_m.a[(burn.in+i),])
zz <- matrix(0,q+1,1)
for(j in 1:n){
if(homo==0){
dfg <- diag(q)
diag(dfg) <- heter$sigma2xi[j,]
sigma_mp[1,1] <- heter$sigma2y[j] + crossprod(rho.a[(burn.in+i),],sigma2_m.a2[,,(burn.in+i)]%*%rho.a[(burn.in+i),])
sigma_mp[2:(q+1),2:(q+1)] <- dfg + sigma2_m.a2[,,(burn.in+i)]
inv.sigma_mp <- solve(sigma_mp)
det_sigma_mp <- det(sigma_mp)
}
zz[1] <- y[j]- comp_s[j]
zz[2:(q+1)] <- M[j,] - mu_m.a[(burn.in+i),]
ss <- sqrt(crossprod(zz,inv.sigma_mp)%*%zz)
D_bar <- D_bar - 2*log(pdf(ss,nu.a[(burn.in+i),])) + log(det_sigma_mp)
mi[j,i] <- det_sigma_mp^(1/2)/pdf(ss,nu.a[(burn.in+i),])
}
D_bar <- D_bar
}
D_bar <- D_bar/post.sam.s


D_theta <- 0
if(q >1){
comp_s <- X%*%apply(as.matrix(beta.a[size,]),2,mean) + B%*%apply(as.matrix(alpha.a[size,]),2,mean) + sum(apply(as.matrix(rho.a[size,]),2,mean)*apply(as.matrix(mu_m.a[(size),]),2,mean))
sigma_mp[1,1] <-  crossprod(apply(as.matrix(rho.a[size,]),2,mean),apply(sigma2_m.a2[,,size],1:2,mean))%*%apply(as.matrix(rho.a[size,]),2,mean)
sigma_mp[1,2:(q+1)] <- crossprod(apply(as.matrix(rho.a[size,]),2,mean),apply(sigma2_m.a2[,,size],1:2,mean))
sigma_mp[2:(q+1),1] <- t(sigma_mp[1,2:(q+1)])
sigma_mp[2:(q+1),2:(q+1)] <- apply(sigma2_m.a2[,,size],1:2,mean)
}else{
comp_s <- X%*%apply(as.matrix(beta.a[size,]),2,mean) + B%*%apply(as.matrix(alpha.a[size,]),2,mean) + sum(apply(as.matrix(rho.a[size,]),2,mean)*apply(as.matrix(mu_m.a[(size),]),2,mean))
sigma_mp[1,1] <- crossprod(apply(as.matrix(rho.a[size,]),2,mean),mean(as.vector(sigma2_m.a2[,,size])))%*%apply(as.matrix(rho.a[size,]),2,mean)
sigma_mp[1,2:(q+1)] <- crossprod(apply(as.matrix(rho.a[size,]),2,mean),mean(as.vector(sigma2_m.a2[,,size])))
sigma_mp[2:(q+1),1] <- t(sigma_mp[1,2:(q+1)])
sigma_mp[2:(q+1),2:(q+1)] <- mean(as.vector(sigma2_m.a2[,,size]))
}
if(homo==1){
sigma_mp[1,1] <- sigma_mp[1,1] + mean(sigma2_y.a[size])
sigma_mp[2:(q+1),2:(q+1)] <- sigma_mp[2:(q+1),2:(q+1)]  + omeg*mean(sigma2_y.a[size])
inv.sigma_mp <- solve(sigma_mp)
det_sigma_mp <- det(sigma_mp)
}
zz <- matrix(0,q+1,1)
mu_M.mean <- apply(as.matrix(mu_m.a[size,]),2,mean)
nu.mean <- apply(as.matrix(nu.a[size,]),2,mean)
for(j in 1:n){
if(homo==0){
dfg <- diag(q)
diag(dfg) <- heter$sigma2xi[j,]
sigma_mp[1,1] <- sigma_mp[1,1] + heter$sigma2y[j]
sigma_mp[2:(q+1),2:(q+1)] <- sigma_mp[2:(q+1),2:(q+1)]  + dfg
inv.sigma_mp <- solve(sigma_mp)
det_sigma_mp <- det(sigma_mp)
}
zz[1] <- y[j]- comp_s[j]
zz[2:(q+1)] <- M[j,] - mu_M.mean
ss <- sqrt(crossprod(zz,inv.sigma_mp)%*%zz)
D_theta <- D_theta - 2*log(pdf(ss,nu.mean)) + log(det_sigma_mp)
}

D_theta <- D_theta 

if(homo==1) res_q <- cdf((y-comp_s)/sqrt(sigma_mp[1,1]), nu.mean)
else res_q <- cdf((y-comp_s)/sqrt(heter$sigma2y), nu.mean)

res_q <- ifelse(res_q < 1e-15, 1e-15, res_q)
res_q <- ifelse(res_q > (1-1e-15), (1-1e-15), res_q)
res_q <- qnorm(res_q)

DIC <- 2*D_bar - D_theta
LMPL <- -sum(log(apply(mi,1,mean)))

KL <- log(apply(mi,1,mean))+ apply(-log(mi),1,mean) #distancia K-L
X_2 <- apply(mi^2,1,mean)/(apply(mi,1,mean))^2- 1  #distancia X^2


list(chains=aa2, DIC=DIC, LMPL=LMPL, residuos=res_q, KL=KL, X_2=X_2)
}
