mcmc.gesm <-
function(params){

rmvnorm.l <- function(mean, sigma){
d <- length(mean)
sigma2 <- chol(sigma)
t(mean + t(sigma2)%*%rnorm(d))
}

MH.gamma <- function(gamma){
y_til <- ((f-1)/4)* Z%*%gamma + (ss*vs - 1)/2

G_b <- solve(S.gamma.II + ((f-1)/4) *ZZ)

g_b <- G_b%*% (S.gamma.II %*%gamma0 + t(Z)%*%y_til)

if(sum(is.na(g_b))>0) g_b <- gamma

d <- length(g_b)
    sigma2 <- chol(G_b)

cont <- 0
    acep <- 0
while(acep==0){
    gamma_new <- g_b + crossprod(sigma2,rnorm(d))
delta <- exp(-0.5*( sum(s2*(exp(-Z%*%gamma_new )-exp(-Z%*%gamma))) + sum(Z%*%(gamma_new - gamma)) + (t(gamma_new - gamma)%*%S.gamma.II %*%(gamma_new - gamma))))
u <- runif(1)
acep <- ifelse(u<delta, 1, 0)
cont <- cont + 1
if(cont == 5){
  ww <- function(gamma_new)(0.5*(sum(s2*(exp(-Z%*%gamma_new ))) + sum(Z%*%(gamma_new)) + t(gamma_new-gamma0)%*%S.gamma.II %*%(gamma_new-gamma0)))
  gamma_new <- optim(gamma,ww,method="BFGS")$par
  #acep <- 1
  g_b <- gamma_new
}
}
gamma_new
}



MH.gamma2 <- function(lambda,gamma){
y_til <- ((f-1)/4)* (D%*%lambda + Z%*%gamma) + (ss*vs - 1)/2

G_b <- solve(prp + ((f-1)/4)*DZDZ)

g_b <- G_b%*%(t(DZ)%*%y_til)

if(sum(is.na(g_b))>0) g_b <- c(gamma,lambda)

d <- length(g_b)
    sigma2 <- chol(G_b)

cont <- 0
    acep <- 0

if(q > 0) gl <- c(gamma,lambda)
else gl <- lambda

while(acep==0){
    gl_new <- g_b + crossprod(sigma2,rnorm(d))
delta <- exp(-0.5*( sum(s2*(exp(-DZ%*%gl_new)-exp(-DZ%*%gl))) + sum(DZ%*%(gl_new - gl)) + (t(gl_new - gl)%*%prp%*%(gl_new - gl))))
u <- runif(1)
acep <- ifelse(u<delta, 1, 0)
cont <- cont + 1
if(cont == 5){
  ww <- function(gl_new)(0.5*(sum(s2*(exp(-DZ%*%gl_new))) + sum(DZ%*%gl_new) + (t(gl_new)%*%prp%*%(gl_new))))
  gl_new <- optim(gl,ww,method="BFGS")$par
  g_b <- gl_new
}
}
if(q >0 ) list(gamma=gl_new[1:q], lambda=gl_new[(1+q):(q+sum(ks2))])
else list(gamma=0, lambda=gl_new[(1+q):(q+sum(ks2))])
}


burn.in <- params$burn.in
post.sam.s <- params$post.sam.s
thin <- params$thin
homo <- params$homo
y <- params$y
p <- params$p
q <- params$q
ks <- params$ks
ks2 <- params$ks2
u <- params$u
pdf <- params$pdf
cdf <- params$cdf
n <- params$n
v <- params$v
fg <- params$fg
kappa <- params$kappa
family <- params$family
total <- burn.in + post.sam.s*thin
ancho <- floor(seq(2, total, length=10))

bar <- txtProgressBar(min=0, max=ancho[10], initial=0, width=50, char="+", style=3)

if(p > 0){
beta.a <- matrix(0,total,p)
beta0 <- matrix(0,p,1)
S.beta <- diag(p)*1000
S.beta.II <- solve(S.beta)
X <- params$X
beta.a[1,] <- params$beta.i}
else{beta.a <- matrix(0,total,1)
     X <- matrix(1,n,1)}

if(sum(ks) > 0){
   alpha.a <- matrix(0,total,sum(ks))
   alpha.a[1,] <- params$alpha.i
   B <- params$B
   tau.a <- matrix(0,total,length(ks))
   a.tau <- 0.001
   b.tau <- 0.001}
else{alpha.a <- matrix(0,total,1)
     tau.a <- matrix(0,total,1)
     B <- matrix(1,n,1)}

if(q > 0){
gamma.a <- matrix(0,total,q)
gamma0 <- matrix(0,q,1)
S.gamma <- diag(q)*1000
S.gamma.II <- solve(S.gamma)
Z <- params$Z
ZZ <- crossprod(Z)
gamma.a[1,] <- params$gamma.i}
else{gamma.a <- matrix(0,total,1)
    Z <- matrix(1,n,1)}


   
if(sum(ks2) > 0){
   lambda.a <- matrix(0,total,sum(ks2))
   lambda.a[1,] <- params$lambda.i
   lambda0 <- matrix(0,sum(ks2),1)
   D <- params$D
   if(q > 0)  DZ <- cbind(Z,D)
   else DZ <- D
   DZDZ <- crossprod(DZ)
   
   tau.l <- matrix(0,total,length(ks2))
   a.tau.l <- 0.001
   b.tau.l <- 0.001}
else{lambda.a <- matrix(0,total,1)
     tau.l <- matrix(0,total,1)
     D <- matrix(1,n,1)}
 
cont <- 1
uii <- matrix(0,n,1)

if(family!="Normal" && family!="Laplace" && attr(params$eta,"know")==0){
extra.parameter <- params$extra.parameter
nu0 <- params$nu0
nu.a <- matrix(0,total,length(nu0))
    nu.a[1,] <- nu0
}else{
nu.a <- kronecker(matrix(1,total,length(params$eta)),t(params$eta))
fknown <- fg(params$eta)
}



for(l in 2:total){

sigma2_y.a <- exp(Z%*%gamma.a[l-1,] + D%*%lambda.a[l-1,])
ss0 <- (y - X%*%beta.a[l-1,] - B%*%alpha.a[l-1,] )/sqrt(sigma2_y.a)
ss <- ss0^2

if(family!="Hyperbolic" && family!="Laplace") uii <- u(ss,nu.a[l-1,])
else for(i in 1:n){uii[i] <- u(ss[i],nu.a[l-1,])}
 
if(family!="Normal" && family!="Laplace"  && attr(params$eta,"know")==0){
nu.a[l,] <- extra.parameter(nu.a[l-1,], uii, ss) 
}

if(p > 0){
sig.beta.a <- solve(S.beta.II + crossprod(X,matrix(1/(sigma2_y.a * kappa(uii)),n,ncol(X))*X))
mu.beta.a <-  sig.beta.a%*%(crossprod(X,matrix(1/(sigma2_y.a*kappa(uii)),n,1)*(y-B%*%alpha.a[l-1,])))
beta.a[l,] <- rmvnorm.l(mean=mu.beta.a, sigma=sig.beta.a)
}
if(sum(ks) > 0){
alpha.a[l,] <- alpha.a[l-1,]

if(length(ks)==1){
scale.tau <- 2*b.tau + crossprod(alpha.a[l-1,],alpha.a[l-1,])
tau.a[l,1] <- 1/rgamma(1, shape= (ks/2 + a.tau), scale= 2/scale.tau)

sig.alpha.a <- solve(diag(ks)/tau.a[l,1] + crossprod(B,matrix(1/(sigma2_y.a*kappa(uii)),n,ks)*B))
mu.alpha.a <- sig.alpha.a%*%(crossprod(B,(1/(sigma2_y.a*kappa(uii)))*(y-X%*%beta.a[l,])))
alpha.a[l,] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)
}
else{
lmi <- 1
for(i in 1:length(ks)){
    B.i <- B[,lmi:(lmi+ks[i]-1)] 
scale.tau <- 2*b.tau + crossprod(alpha.a[l-1,lmi:(lmi+ks[i]-1)],alpha.a[l-1,lmi:(lmi+ks[i]-1)])
tau.a[l,i] <- 1/rgamma(1, shape= (ks[i]/2 + a.tau), scale= 2/scale.tau)

sig.alpha.a <- solve(diag(ks[i])/tau.a[l,i] + crossprod(B.i,matrix(1/(sigma2_y.a*kappa(uii)),n,ks[i])*B.i))
mu.alpha.a <- sig.alpha.a%*%(crossprod(B.i,(1/(sigma2_y.a*kappa(uii)))*(y-X%*%beta.a[l,] - B[,-c(lmi:(lmi+ks[i]-1))]%*%alpha.a[l,-c(lmi:(lmi+ks[i]-1))])))
alpha.a[l,lmi:(lmi+ks[i]-1)] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)

lmi <- lmi + ks[i]
}
}
}

s2 <- (y - X%*%beta.a[l,] - B%*%alpha.a[l,] )^2/kappa(uii)
ss <- (y - X%*%beta.a[l,] - B%*%alpha.a[l,] )^2/sigma2_y.a

if(homo==1){
scale.sigma2_y <- sum(s2) + 0.0001/2
gamma.a[l] <- log(1/rgamma(1,shape= ((n + 0.0001/2)/2), scale=2/scale.sigma2_y))
}else{
vs <- v(sqrt(ss),nu.a[l,])
if(family=="Normal" || family=="Laplace" || attr(params$eta,"know")==1) f <- fknown
else f <- fg(nu.a[l,])

if(sum(ks2) >0){
prp <- matrix(0,q+sum(ks2),q+sum(ks2))
lim <- cumsum(c(0,ks2))
if(q > 0) prp[1:q,1:q] <- S.gamma.II

for(i in 1:length(ks2)) {
scale.tau.l <- 2*b.tau.l + sum((lambda.a[l-1,(lim[i]+1):lim[i+1]])^2)
tau.l[l,i] <- 1/rgamma(1, shape= (ks2[i]/2 + a.tau.l), scale= 2/scale.tau.l)
prp[(lim[i]+1+q):(lim[i+1]+q),(lim[i]+1+q):(lim[i+1]+q)] <- diag(ks2[i])/tau.l[l,i]
}
MH <- MH.gamma2(lambda.a[l-1,],gamma.a[l-1,])
gamma.a[l,] <- MH$gamma
lambda.a[l,] <- MH$lambda
}
else{
gamma.a[l,] <- MH.gamma(gamma.a[l-1,])
    sigma2_y.a <- exp(Z%*%gamma.a[l,] + D%*%lambda.a[l-1,])
}
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
aa <- matrix(0,post.sam.s,p+q+sum(ks)+sum(ks2)+length(nu0)+length(ks)+length(ks2))}
else{aa <- matrix(0,post.sam.s,p+q+sum(ks)+sum(ks2)+length(ks)+length(ks2))}

cad <- as.vector(" ")

if(p > 0){
aa[,1:p] <- beta.a[size,]
cad <- cbind(cad,t(paste("beta",1:p)))
}

if(sum(ks) > 0){
aa[,(p+1):(p+sum(ks))] <- alpha.a[size,]
aa[,(p+sum(ks)+1):(p+sum(ks)+length(ks))] <- tau.a[size,]
cad <- cbind(cad,t(colnames(B)))
cad <- cbind(cad,t(paste("tau_alpha",1:length(ks))))
}


if(q > 0){
aa[,(p+sum(ks)+length(ks)+1):(p+sum(ks)+length(ks)+q)] <- gamma.a[size,]
cad <- cbind(cad,t(paste("gamma",1:q)))
}

if(sum(ks2) > 0){
aa[,(p+sum(ks)+length(ks)+q+1):(p+sum(ks)+length(ks)+q+sum(ks2))] <- lambda.a[size,]
aa[,(p+sum(ks)+length(ks)+q+sum(ks2)+1):(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2))] <- tau.l[size,]
cad <- cbind(cad,t(colnames(D)))
cad <- cbind(cad,t(paste("tau_lambda",1:length(ks2))))
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
aa[,(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+1):(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+length(nu0))] <- nu.a[size,]
}


if(family!="Normal" && family!="Laplace"  && attr(params$eta,"know")==0){
cad <- cbind(cad,t(paste("eta",1:length(nu0))))
}

haver <- apply(aa,2,var)
aa2 <- aa[,haver!=0]
cad <- cad[2:length(cad)]
colnames(aa2) <- cad


mi <- matrix(0,n,post.sam.s)
res_q <- matrix(0,n,post.sam.s)

D_bar <- 0
for(i in 1:post.sam.s){
mui <- (X%*%beta.a[(burn.in+i),]+B%*%alpha.a[(burn.in+i),])
sigma2i <- exp(Z%*%gamma.a[(burn.in+i),] + D%*%lambda.a[(burn.in+i),])
zii <- (y-mui)/sqrt(sigma2i)
D_bar <- D_bar - 2*sum(log(pdf(zii,nu.a[(burn.in+i),])) -log(sigma2i)/2)
    mi[,i] <- 1/(pdf(zii,nu.a[(burn.in+i),])/sqrt(sigma2i))
res_q[,i] <- qnorm(cdf(zii,nu.a[(burn.in+i),]))
}
D_bar <- D_bar/post.sam.s


D_theta <- 0
mu_bar <- X%*%apply(as.matrix(beta.a[size,]),2,mean) + B%*%apply(as.matrix(alpha.a[size,]),2,mean)
sigma2_bar <- exp(Z%*%apply(as.matrix(gamma.a[size,]),2,mean) + D%*%apply(as.matrix(lambda.a[size,]),2,mean))
zii <- (y-mu_bar)/sqrt(sigma2_bar)
D_theta <- -2*sum(log(pdf(zii,nu.a[(burn.in+i),])) -log(sigma2_bar)/2)

DIC <- 2*D_bar - D_theta
LMPL <- -sum(log(apply(mi,1,mean)))
res <- apply(res_q,1,mean)            #residuos
KL <- log(apply(mi,1,mean))+ apply(-log(mi),1,mean) #distancia K-L
X_2 <- apply(mi^2,1,mean)/(apply(mi,1,mean))^2- 1  #distancia X^2


list(DIC=DIC, LMPL=LMPL, chains=aa2, residuos=res, KL=KL, X_2=X_2)
}
