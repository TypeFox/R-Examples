prm <-
function(X,y,a,fairct=4,opt="l1m",usesvd=FALSE){

# PRM Partial Robust M-regression estimator
# This version does not rely on mvr{chemometrics} which is always 
# mean-centering X; instead, "Unisimpls" is used, see below

##################################################################
Unisimpls <- function(X,y,a){
  # univariate simpls; PF, 09.02.2010
  # X ... explanatory variables
  # y ... response variable
  # a ... number of desired components
  ##################################
  n <-  nrow(X)
  px <-  ncol(X)
  # if n<px do SVD:
  if (px > n) {
    dimensions <- 1
    dimension <- px - n
    ressvd <- svd(t(X))
    X <- ressvd$v %*% diag(ressvd$d)
    n <- nrow(X)
    px <- ncol(X)
  }
  else {
    dimensions <- 0
  }
  # end of SVD
  s <- t(X)%*%y
  U <- matrix(0,nrow=n,ncol=a)
  V <- matrix(0,nrow=px,ncol=a)
  R <- matrix(0,nrow=px,ncol=a) # simpls weights
  B <- V
  for (j in 1:a) {  # extract one component after the other
    r <- s
    u <- X%*%r
    u <- u-U[,1:max(1,j-1)]%*%(t(U[,1:max(1,j-1)])%*%u)
    normu <- drop(sqrt(t(u)%*%u))
    u <- u/normu
    r <- r/normu
    p <- t(X)%*%u
    v <- p - V[,1:max(1,j-1)]%*%(t(V[,1:max(1,j-1)])%*%p)
    v <- v/drop(sqrt(t(v)%*%v))
    s <- s-v%*%(t(v)%*%s)
    U[,j] <- u
    R[,j] <- r
    V[,j] <- v
    B[,j] <- R[,1:j]%*%t(R[,1:j])%*%t(X)%*%y
  }
  # if SVD:
  if (dimensions == 1) {
    B <- ressvd$u %*% B
  }
  # end of SVD
  list(coefficients=B,scores=U)
}
##################################################################

n=nrow(X)
p=ncol(X)

if (usesvd==TRUE){
  if (p>n){
	dimensions <- 1
	dimension <- p-n
	ressvd <- svd(t(X))
	X <- ressvd$v%*%diag(ressvd$d)
	n=nrow(X)
	p=ncol(X)
  }
  else {dimensions <- 0}
}

if (opt=="l1m"){
	#require(pcaPP)
	mx <- l1median(X)   
}
else if (opt=="mean"){
	mx <- apply(X,2,mean)   
}
else { mx <- apply(X,2,median)}

if (opt=="mean"){
  my <- mean(y)
}
else {
  my <- median(y)
}

Xmc <- as.matrix(scale(X,center=mx,scale=FALSE))
ymc <- as.vector(scale(y, center = my, scale = FALSE))
wx <- sqrt(apply(Xmc^2,1,sum))
wx <- wx/median(wx)
wx <- 1/((1+abs(wx/fairct))^2)

wy <- abs(ymc)
wy <- wy/median(wy)
wy <- 1/((1+abs(wy/fairct))^2)

w <- wx*wy
Xw <- Xmc*sqrt(w)
yw <- ymc*sqrt(w)

loops <- 1
ngamma <- 10^5
difference <- 1

#require(pls)
while ((difference>0.01) && loops<30){
	ngammaold <- ngamma
	###spls <- mvr(yw~Xw,ncomp=a,method="simpls")
        spls <- Unisimpls(Xw,yw,a) # new
	###b <- spls$coef[,,a]
        b <- spls$coef[,a] # new
	###gamma <- t(t(ymc)%*%spls$sco) # wrong
        gamma <- t(t(yw) %*% spls$sco)
	T <- spls$sco/sqrt(w)
	r <- ymc-T%*%gamma
	rc <- r-median(r)
	r <- rc/median(abs(rc))
	wy <- 1/((1+abs(r/fairct))^2)
	if (opt=="l1m"){mt <- l1median(T)}
	else {mt <- apply(T,2,median)}
	dt <- T-mt
	wt <- sqrt(apply(dt^2,1,sum))
	wt <- wt/median(wt)
	wt <- 1/((1+abs(wt/fairct))^2)
	ngamma <- sqrt(sum(gamma^2))
	difference <- abs(ngamma-ngammaold)/ngamma
	w <- drop(wy*wt)
        # new BL: 17.9.09
        w0 <- which(w==0)
        if(length(w0) != 0) { w <- replace(w, list=w0, values=10^(-6)) }
        # corrected:
         Xw <- Xmc*sqrt(w)
         yw <- ymc*sqrt(w)
        #1#Xw <- X * sqrt(w)
        #1#yw <- y * sqrt(w)
	#print(difference)
	#print(loops)
	loops <- loops+1
	}

if (usesvd==TRUE){
  #b0 <- drop(coef(spls, intercept=TRUE))[1]
  #b0 <- drop(t(mx)%*%b)
  if (dimensions==1){
	b <- drop(ressvd$u%*%b)
        if (opt=="mean"){
          b0 <- mean(y-as.matrix(X)%*%t(ressvd$u)%*%b)
        }
        else {
          b0 <- median(y-as.matrix(X)%*%t(ressvd$u)%*%b)
        }
	yfit <- as.matrix(X)%*%t(ressvd$u)%*%b+b0
  }
  else {
	#yfit <- as.matrix(Xmc)%*%b+my+b0
        if (opt=="mean"){
          b0 <- mean(y-as.matrix(X)%*%b)
        }
        else {
          b0 <- median(y-as.matrix(X)%*%b)
        }
	yfit <- as.matrix(X)%*%b+b0
  }
}
else {
        #b0 <- drop(coef(spls, intercept=TRUE))[1]
        ###b0 <- drop(t(mx)%*%b)+my
        ###b0 <- -drop(t(mx)%*%b)+my
        ###yfit <- as.matrix(Xmc)%*%b+b0
        ###yfit <- as.matrix(X)%*%b+b0
        if (opt=="mean"){
          b0 <- mean(y-as.matrix(X)%*%b)
        }
        else {
          b0 <- median(y-as.matrix(X)%*%b)
        }
        yfit <- as.matrix(X)%*%b+b0
        
}

list(coef=b,intercept=b0,wy=wy,wt=wt,w=w,scores=T,loadings=spls$loadings,
	fitted.values=yfit,mx=mx,my=my)
}

