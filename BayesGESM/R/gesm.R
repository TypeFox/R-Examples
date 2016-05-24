gesm <-
function(formula, data, family, eta, burn.in, post.sam.s, thin){

if(family!="Normal" & family!="Student-t" & family!="Slash" & family!="Hyperbolic" & family!="ContNormal" & family!="Laplace")
stop("Family of Distributions specified by the user is not supported, Check documentation!!",call.=FALSE)


if(family=="Laplace" | family=="Normal"){
 if(!missingArg(eta)) stop("for the Laplace and Normal distribution must be not specify the extra parameter!!", call.=FALSE)
 eta <- 0
 attr(eta,"know") <- 1
 }

if(missingArg(eta)){
  if(family=="ContNormal"){eta <- c(0,0)}else{eta <- 0}
  attr(eta,"know") <- 0
}else{attr(eta,"know") <- 1}

if(missingArg(thin)){thin <- 1}

if(thin<1 | thin != floor(thin))
stop("the thin value must be a positive integer!!", call.=FALSE)
if (post.sam.s != floor(post.sam.s) | post.sam.s < 1  ){ stop("Invalid posterior sample size value", call.=FALSE)}
if (burn.in != floor(burn.in) | burn.in < 1){ stop("Invalid burn-in value", call.=FALSE)}


 mf <- match.call(expand.dots = FALSE)
 m <- match(c("formula"), names(mf), 0)
 mf <- mf[c(1, m)]

 if (missingArg(data)) 
 data <- environment(formula)

 formula <- Formula(formula)
 if (length(formula)[2L] < 2L)
        formula <- as.Formula(formula(formula), ~1)

nl <- formula(formula, lhs=0, rhs=1)
ll <- formula(formula, lhs=1, rhs=2)

response <- as.matrix(eval(ll[[2]], data))
n <- length(response)
if(ncol(as.matrix(response)) > 1) stop("The response variable must be univariate!!", call.=FALSE)

 par_ini <- list(family=family,n=n,y=response)
 x <- model.matrix(nl, data)
 xa <- "bsp("
 xb <- colnames(x)
 idx <- grepl(xa,xb,fixed=TRUE)
 
 if(sum(idx) >= 1){ 
    ks <- matrix(0,sum(idx),1)
bss <- matrix(0,nrow(x),1)
cont <- 1
alphas <- " "
for(i in 1:length(idx)){
if(idx[i]==1){
temp <- eval(parse(text=xb[i]), data)
  ks[cont] <- attr(temp,"kn") + 3
        bss <- cbind(bss,attr(temp,"B"))
alphas <- cbind(alphas,t(paste("alpha",cont,1:ks[cont])))
            cont <- cont + 1
}
}
bss <- bss[,-1]
     colnames(bss) <- alphas[-1]
par_ini$B <- bss

nps <- as.matrix(x[,idx==1])
colnames(nps) <- colnames(x)[idx==1]
par_ini$nps <- nps
     p <- sum(!idx)
if(p > 0){
   X <- as.matrix(x[,!idx])
    colnames(X) <- colnames(x)[!idx]
X_au <- cbind(X,bss)
            b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
        par_ini$beta.i <- b_au[1:p]
        par_ini$alpha.i <- b_au[(p+1):(p+sum(ks))]
par_ini$X <- X
}else{X_au <- bss
      b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
      par_ini$alpha.i <- b_au[1:length(ks)]
}
 }else{
    ks <- 0
    X <- as.matrix(x[,])
    colnames(X) <- colnames(x)
    p <- ncol(X)
X_au <- X
        b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
    par_ini$beta.i <- b_au
par_ini$X <- X
 }
 par_ini$p <- p
 par_ini$ks <- ks
 rres <- log((response -  X_au%*%b_au)^2)


 homo <- 0 
 z <- model.matrix(ll, data)
 za <- "bsp("
 zb <- colnames(z)
 idz <- grepl(za,zb,fixed=TRUE)
 
 if(sum(idz) >= 1){ 
    ks2 <- matrix(0,sum(idz),1)
bss <- matrix(0,nrow(z),1)
cont <- 1
lambdas <- " "
for(i in 1:length(idz)){
if(idz[i]==1){
temp <- eval(parse(text=zb[i]), data)
  ks2[cont] <- attr(temp,"kn") + 3
        bss <- cbind(bss,attr(temp,"B"))
lambdas <- cbind(lambdas,t(paste("lambda",cont,1:ks2[cont])))
            cont <- cont + 1
}
}
bss <- bss[,-1]
colnames(bss) <- lambdas[-1]
par_ini$D <- bss

nps2 <- as.matrix(z[,idz==1])
colnames(nps2) <- colnames(z)[idz==1]
par_ini$nps2 <- nps2

     q <- sum(!idz)
if(q > 0){
   Z <- as.matrix(z[,!idz])
    colnames(Z) <- colnames(z)[!idz]
Z_au <- cbind(Z,bss)
            b_au <- solve(t(Z_au)%*%Z_au)%*%t(Z_au)%*%rres
        par_ini$gamma.i <- b_au[1:q]
        par_ini$lambda.i <- b_au[(q+1):(q+sum(ks2))]
par_ini$Z <- Z
}else{Z_au <- bss
      b_au <- solve(t(Z_au)%*%Z_au)%*%t(Z_au)%*%rres
      par_ini$lambda.i <- b_au[1:length(ks2)]
}
 }else{
    ks2 <- 0
    Z <- as.matrix(z[,])
    colnames(Z) <- colnames(z)
    q <- ncol(Z)
if(q==1 && sd(Z[,1])==0) homo <- 1
Z_au <- Z
        b_au <- solve(t(Z_au)%*%Z_au)%*%t(Z_au)%*%rres
    par_ini$gamma.i <- b_au
par_ini$Z <- Z
 }

 par_ini$q <- q
 par_ini$ks2 <- ks2

 if(family=="Normal"){
    u <- function(s,eta) matrix(1,length(s),1) 
v <- function(z,eta) matrix(1,length(z),1)
fg <- function(eta) 3

pdf <- function(z,eta){dnorm(z)}
cdf <- function(z,eta){pnorm(z)}

}
 
 if(family=="Student-t"){
 if(attr(eta,"know") == 1){
    if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)}
 
    u <- function(s,eta){
    bb <- eta/2 + s/2
 aa <- (eta + 1)/2 
rgamma(length(s),shape=aa,scale=1)/bb
    }
v <- function(z,eta) (eta+1)/(eta + z^2)
fg <- function(eta) 3*(eta+1)/(eta+3)

  pdf <- function(z,eta){dt(z,eta)}
 cdf <- function(z,eta){pt(z,eta)}

  extra.parameter <- function(nu0, ui,s){
       av <- 0.01
       bv <- 0.01
       log.poster <- function(nu){((length(ui)*nu + 2*av - 2)/2)*log(nu/2) - length(ui)*log(gamma(nu/2)) - nu*(sum(ui)-sum(log(ui)) + 2*bv)/2}
       tf <- function(a){
       nu <- exp(a)
       -log.poster(nu)
   }
       out <- optim(log(nu0), tf, method="L-BFGS-B", lower=log(1), upper=log(40))
       new <- exp(out$par)
       new
}
     nu0 <- 5
 }

 
 if(family=="Slash"){
  if(attr(eta,"know") == 1){
 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)}
 
     u <- function(s,eta){
     bb <- s/2
     aa <- 1/2 + eta 
 u <- runif(length(s))*pgamma(bb, shape=aa, scale=1)
 ui <- qgamma(u, shape=aa, scale=1)/bb
 ifelse(ui< 1e-20,1e-20,ui)
   }

   
G <- function(a,x) gamma(a)*pgamma(1,shape=a,scale=1/x)/(x^a)
v <- function(z,eta) G(eta+3/2,z^2/2)/G(eta+1/2,z^2/2)   

ds <- function(z,eta) eta*G(eta+1/2,z^2/2)/sqrt(2*pi)
gfg <- function(z,eta) ds(z,eta)*(v(z,eta))^2*z^4
fg <- function(eta) 2*integrate(gfg,0,Inf,eta)$value
  

       pdf <- function(z,eta){eta*(z^2/2)^(-(eta+1/2))*gamma(eta+1/2)*pgamma(1,shape=(eta+1/2), scale=(2/z^2))/(2*pi)^(1/2)}
    cdf <- function(z,eta){
   temp <- matrix(0,length(z),1)          
   for(i in 1:length(z)) temp[i] <- integrate(pdf,-Inf,z[i],eta)$value    
   temp    
}   

  extra.parameter <- function(nu0, ui, s){
       av <- 0.1
       bv <- 0.1

   new <- rgamma(1,shape=(length(ui)+av), scale= (1/(bv-sum(log(ui)))))
   min(new,40)
       }
   nu0 <- 8
   }
 
 if(family=="Laplace"){
     u <- function(s,eta){
     bb <- s
 aa <- 1/2
     rgig(1,lambda=aa,chi=bb,psi=0.25)
 }

 v <- function(z,eta) abs(z)^(-1/2)
 fg <- function(eta) 2
 
  pdf <- function(z,eta){dnormp(z,mu=0,sigmap=2,p=1)}
 cdf <- function(z,eta){pnormp(z,mu=0,sigmap=2,p=1)}

 }

 if(family=="Hyperbolic"){
 if(attr(eta,"know") == 1){
   if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)}
  
    u <- function(s,eta){
    bb <- s + 1
aa <-  1/2
      rgig(1,lambda=aa,chi=bb,psi=eta^2)
    }

 v <- function(z,eta) eta/sqrt(1 + z^2)
 dh <- function(z,eta) exp(-eta*sqrt(1+z^2))/(2*besselK(eta, 1))

 fgf <- function(z,eta) dh(z,eta)*(eta*z^2/sqrt(1 + z^2))^2
     fg <- function(eta) 2*integrate(fgf,0,Inf,eta)$value

pdf <- function(z,eta) dh(z,eta)
cdf <- function(z,eta){
temp <- matrix(0,length(z),1)
for(i in 1:length(z)) temp[i] <- integrate(dh,-Inf,z[i],eta)$value
temp/(2*integrate(dh,0,Inf,eta)$value)
}


extra.parameter <- function(nu0, ui, s){
       av <- 0.0001
       bv <- 0.0001
       log.poster <- function(nu){(length(ui)+av - 1)*log(nu) - length(ui)*log(besselK(nu,1)) - (nu^2*sum(ui) + 2*bv * nu)/2}

       tf <- function(a){
        nu <- exp(a)
       -log.poster(nu)
   }
       out <- optim(log(nu0), tf, method="L-BFGS-B", lower=log(0.2), upper=log(10))
       new <- exp(out$par)
       new
}
     nu0 <- 1
 }

 
 if(family=="ContNormal"){
 if(attr(eta,"know") == 1){
 if(eta[1]<0 | eta[1]>1) stop("the extra parameter eta[1] must be within the interval (0,1)!!",call.=FALSE)
 if(eta[2]<0 | eta[2]>1) stop("the extra parameter eta[2] must be within the interval (0,1)!!",call.=FALSE)}
 
 u <- function(s,eta){
 bb <- exp(-s*eta[2]/2)*eta[1]*eta[2]^(1/2)
 aa <- (1 - eta[1])*exp(-s/2)
 bb <- ifelse(aa==0 && bb==0, 1, bb/(aa + bb))
 uu <- runif(length(s))
 uii <- ifelse(uu<=bb,eta[2],1)
 uii
 }

 v <- function(z,eta) (eta[2]^(3/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))/(eta[2]^(1/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))
    fgd <- function(z) (v(z,eta))^2*z^4*(sqrt(eta[2])*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))
fg <- function(eta) 2*integrate(fgd,0,35)$value


 pdf <- function(z,eta){eta[1]*dnorm(z,sd=1/sqrt(eta[2]))+(1 - eta[1])*dnorm(z)}
 cdf <- function(z,eta){eta[1]*pnorm(z,sd=1/sqrt(eta[2]))+(1 - eta[1])*pnorm(z)}
 
   extra.parameter <- function(nu0, ui, s){
   av1 <- 0.0001
   bv1 <- 0.0001
   av2 <- 2
   bv2 <- 2
   uis <- ifelse(ui==1,0,1)
   nv2 <- sum(uis)

      aa <- (nv2*+av2)/2
   bb <- (sum(s*uis/2) + bv2)
   u <- runif(1)*pgamma(1, shape=aa, scale=1/bb)

   new1 <- max(min(rbeta(1, shape1=(nv2+av1), shape2=(length(ui)-nv2+bv1)),0.99),0.01)
   new2 <- min(max(qgamma(u, shape=aa, scale=1/bb),0.01),0.99)
   c(new1,new2)
}
nu0 <- c(0.5,0.5)
 }
 
 if(family=="Student-t" | family=="Slash" | family=="ContNormal"){
 kappa <- function(u){
 1/u  }
 }
 else{kappa <- function(u){
 u  } }

 if(family!="Normal" && family!="Laplace"){
 par_ini$extra.parameter <- extra.parameter
    par_ini$nu0 <- nu0
 }

   par_ini$u <- u
   par_ini$burn.in <- burn.in
   par_ini$post.sam.s <- post.sam.s
   par_ini$thin <- thin
   par_ini$kappa <- kappa
   par_ini$eta <- eta
   par_ini$pdf <- pdf
   par_ini$cdf <- cdf
   par_ini$homo <- homo
   par_ini$v <- v
   par_ini$fg <- fg
   par_ini$sigma2_y <- rres
   result <- mcmc.gesm(par_ini)
par_ini$chains=result$chains
par_ini$DIC=result$DIC
par_ini$LMPL=result$LMPL
par_ini$res=result$res
par_ini$KL=result$KL
par_ini$X_2=result$X_2
class(par_ini) <- "gesm"
par_ini$call <- match.call()
par_ini
}
