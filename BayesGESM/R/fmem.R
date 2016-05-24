fmem <-
function(formula, data, omeg, family, eta, burn.in, post.sam.s, thin, heter){

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

 x <- model.matrix(nl, data)
 if(length(x)==0) stop("At least one covariate with measurement error must be specified!!", call.=FALSE) 

 M <- as.matrix(x[,-1])
 colnames(M) <- colnames(x)[-1]
 xa <- "bsp("
 xb <- colnames(M)
 q <- ncol(M)
 idx <- grepl(xa,xb,fixed=TRUE)
 if(sum(idx) >= 1) stop("Nonlinear effects of  covariates with measurement error  are not supported!!",call.=FALSE)


 if(missingArg(heter)){
   homo <- 1
   if(missingArg(omeg)){omeg <- 1}
   if(omeg<=0){stop("The value of the Ratio of the error variances must be positive", call.=FALSE)}
}else{
homo <- 0
if(!is.list(heter)) stop("The argument must be a list!!", call.=FALSE)

heter$sigma2y <- as.matrix(heter$sigma2y)
heter$sigma2xi <- as.matrix(heter$sigma2xi)
if(nrow(heter$sigma2y) != n || ncol(heter$sigma2y)!= 1 || nrow(heter$sigma2xi) != n || ncol(heter$sigma2xi)!= q) stop("The matrix for heteroscedastic model are ...")
}


 w <- model.matrix(ll, data)
 wa <- "bsp("
 wb <- colnames(w)
 idw <- grepl(wa,wb,fixed=TRUE)
 
 if(sum(idw) >= 1){ 
    ks <- matrix(0,sum(idw),1)
bss <- matrix(0,nrow(w),1)
cont <- 1
alphas <- " "
for(i in 1:length(idw)){
if(idw[i]==1){
temp <- eval(parse(text=wb[i]), data)
  ks[cont] <- attr(temp,"kn") + 3
        bss <- cbind(bss,attr(temp,"B"))
alphas <- cbind(alphas,t(paste("alpha",cont,1:ks[cont])))
            cont <- cont + 1
}
}
bss <- bss[,-1]
   X <- as.matrix(w[,!idw])
nps <- as.matrix(w[,idw==1])
    colnames(X) <- colnames(w)[!idw]
    colnames(nps) <- colnames(w)[idw==1]
     p <- sum(!idw)
 }else{
    ks <- 0
    X <- as.matrix(w[,])
    colnames(X) <- colnames(w)
    p <- ncol(X)
 }

par_ini <- list(p=0,ks=ks,q=0,family=family,n=n,y=response)
par_ini$q <- q
par_ini$M <- M

     nombres <- colnames(X)
     if(sum(ks)==0){
 X_au <- cbind(X,M)
 if(homo==0){
 PP <- matrix(1/heter$sigma2y, nrow(X_au), ncol(X_au))
 b_au <- solve(t(X_au)%*%(X_au*PP))%*%t(X_au*PP)%*%response
 }
 else b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
 
 par_ini$p <- p
 par_ini$X <- X
 par_ini$beta.i <- b_au[1:p]
 par_ini$rho.i <- b_au[(p+1):(q+p)]
 
 rres <- mean((response -  X_au%*%b_au)^2)
 }
 else{
 B <- bss
  colnames(B) <- alphas[-1]
  par_ini$p <- p
  par_ini$X <- X
  par_ini$B <- B
  par_ini$nps <- nps
  X_au <- cbind(X,B,M)
  if(homo==0){
 PP <- matrix(1/heter$sigma2y, nrow(X_au), ncol(X_au))
 b_au <- solve(t(X_au)%*%(X_au*PP))%*%t(X_au*PP)%*%response
  } 
  else b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
  par_ini$beta.i <- b_au[1:p]
  par_ini$alpha.i <- b_au[(p+1):(p+sum(ks))]
  par_ini$rho.i <- b_au[(p+sum(ks)+1):(p+sum(ks)+q)]
  rres <- mean((response -  X_au%*%b_au)^2)
      }  
 
 if(family=="Normal"){
    u <- function(s,eta){
    matrix(1,length(s),1) 
    }
 pdf <- function(z,eta){dnorm(z)/(2*pi)^(q/2)}
  cdf <- function(z,eta){pnorm(z)}
 }
 
 if(family=="Student-t"){
 if(attr(eta,"know") == 1){
 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)}
 
    u <- function(s,eta){
    bb <- eta/2 + s/2
 aa <- (eta + 1)/2 + q
rgamma(length(s),shape=aa,scale=1)/bb
    }

 pdf <- function(z,eta){gamma((q+1+eta)/2)*(1+z^2/eta)^(-(q+1+eta)/2)/(gamma(eta/2)*(pi*eta)^((q+1)/2))}
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
     aa <- 1/2 + eta + q
 u <- runif(length(s))*pgamma(bb, shape=aa, scale=1)
 qgamma(u, shape=aa, scale=1)/bb
   }

 pdf <- function(z,eta){eta*(z^2/2)^(-(eta+(q+1)/2))*gamma(eta+(q+1)/2)*pgamma(1,shape=(eta+(q+1)/2), scale=(2/z^2))/(2*pi)^((q+1)/2)}
 pdf2 <- function(z,eta){eta*(z^2/2)^(-(eta+(1)/2))*gamma(eta+(1)/2)*pgamma(1,shape=(eta+(1)/2), scale=(2/z^2))/(2*pi)^((1)/2)}

 cdf <- function(z,eta){temp <- matrix(0,length(z),1)          
 for(i in 1:length(z)){    
 temp[i] <- integrate(pdf2,-Inf,z[i],eta)$value    
 }    
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
 aa <- 1/2 -q
     rgig(1,lambda=aa,chi=bb,psi=0.25)
 }
 pdf <- function(z,eta){besselK(sqrt(z^2/4),(-(q+1)/2+1))*(z^2)^(-((q+1)-2)/4)/((pi)^((q+1)/2))*4^(-(q+2)/2)}

 cdf <- function(z,eta){pnormp(z,mu=0,sigmap=2,p=1)}
 }

 if(family=="Hyperbolic"){
 if(attr(eta,"know") == 1){
   if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)}
  
    u <- function(s,eta){
    bb <- s + 1
aa <- -q + 1/2
      rgig(1,lambda=aa,chi=bb,psi=eta^2)
    }

pdf <- function(z,eta){besselK((sqrt(z^2+1)*eta),(-(q+1)/2+1))*eta^((q+1)/2)*(z^2+1)^((-(q+1)+2)/4)/((2*pi)^((q+1)/2)*besselK(eta,1))}
pdf2 <- function(z,eta){besselK((sqrt(z^2+1)*eta),(-(1)/2+1))*eta^((1)/2)*(z^2+1)^((-(1)+2)/4)/((2*pi)^((1)/2)*besselK(eta,1))}

cdf <- function(z,eta){
temp <- matrix(0,length(z),1)
for(i in 1:length(z)){
            temp[i] <- integrate(pdf2,-Inf,z[i],eta)$value
}
temp/(2*integrate(pdf2,0,Inf,eta)$value)
}
    extra.parameter <- function(nu0, ui, s){
       av <- 0.0001
       bv <- 0.0001
       log.poster <- function(nu){(length(ui)+av - 1)*log(nu) - length(ui)*log(besselK(nu,1)) - (nu^2*sum(ui) + 2*bv * nu)/2}

       tf <- function(a){
        nu <- exp(a)
       -log.poster(nu)
   }
       out <- optim(log(nu0), tf, method="L-BFGS-B", lower=log(0.2), upper=log(20))
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
 bb <- exp(-s*eta[2]/2)*eta[1]*eta[2]^(1/2 + q)
 aa <- (1 - eta[1])*exp(-s/2)
 bb <- ifelse(aa==0 && bb==0, 1, bb/(aa + bb))
 uu <- runif(length(s))
 uii <- ifelse(uu<=bb,eta[2],1)
 uii
 }

 pdf <- function(z,eta){eta[2]^((q+1)/2)*eta[1]*dnorm(z*sqrt(eta[2]))/(2*pi)^(q/2) + (1-eta[1])*dnorm(z)/(2*pi)^(q/2)}  

 cdf <- function(z,eta){eta[1]*pnorm(z,sd=1/sqrt(eta[2]))+(1 - eta[1])*pnorm(z)}

   extra.parameter <- function(nu0, ui, s){
   av1 <- 0.0001
   bv1 <- 0.0001
   av2 <- 2
   bv2 <- 2
   nv2 <- sum(ifelse(ui==1,0,1))

   aa <- (nv2*(q/2+1)+av2)
   bb <- (sum(s*ifelse(ui==1,0,1)/2) + bv2)
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
   par_ini$sigma2_y <- rres
   par_ini$homo <- homo
   if(homo == 0) par_ini$heter <- heter
   else  par_ini$omeg <- omeg

   result <- mcmc.fmem(par_ini)
   par_ini$chains=result$chains
   par_ini$DIC=result$DIC
   par_ini$LMPL=result$LMPL
   par_ini$res=result$residuos
   par_ini$KL=result$KL
   par_ini$X_2=result$X_2
   class(par_ini) <- "fmem"
   par_ini$call <- match.call()
   par_ini
   }
