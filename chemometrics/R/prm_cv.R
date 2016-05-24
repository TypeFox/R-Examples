prm_cv <-
function(X,y,a,fairct=4,opt="median",subset=NULL,
	segments=10,segment.type="random",trim=0.2,sdfact=2,plot.opt=TRUE){

# PRM Partial Robust M-regression estimator with CV
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
} # end of Unisimpls
##################################################################


##################################################################
# calls pcr_cvsub for cross validating prm (partial robust M-regression)
prm_cvsub <- function(X,y,a,fairct,opt,subset,
	segments,segment.type,trim){

if (!is.null(subset)){
	X <- X[subset,]
	y <- y[subset]
}

n=nrow(X)
p=ncol(X)

if (opt=="l1m"){
	#require(pcaPP)
	mx <- l1median(X)   
}
else { mx <- apply(X,2,median)}
my <- median(y)

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

#require(pls)

# start CV
ypred <- rep(NA,length(y))
segment <- cvsegments(n, k = segments, type = segment.type)
SEPj <- rep(NA,segments) # SEP in each segment
SEPtrimj <- rep(NA,segments) # SEP trimmed in each segment
for (n.seg in 1:length(segment)) {
	seg <- segment[[n.seg]]
        obsuse <- as.numeric(unlist(segment[-n.seg]))
	wsel <- w[obsuse]
	Xwsel <- Xw[obsuse,]
	ywsel <- yw[obsuse]
	dsel=list(ywsel=ywsel,Xwsel=as.matrix(Xwsel))
        loops <- 1
        ngamma <- 10^5
        difference <- 1
        while ((difference>0.01) && loops<30){
          ngammaold <- ngamma
          ###spls <- mvr(ywsel~Xwsel,data=dsel,ncomp=a,method="simpls")
	  spls <- Unisimpls(Xwsel,ywsel,a) # new
	  ###b <- spls$coef[,,a]
	  b <- spls$coef[,a] # new
	  ###gamma <- t(t(ymc[obsuse])%*%spls$sco) # wrong
	  gamma <- t(t(ywsel) %*% spls$sco)
	  T <- spls$sco/sqrt(wsel)
	  r <- ymc[obsuse]-T%*%gamma
	  rc <- r-median(r)
	  r <- rc/median(abs(rc))
	  wysel <- 1/((1+abs(r/fairct))^2)
	  if (opt=="l1m"){mt <- l1median(T)}
	  else {mt <- apply(T,2,median)}
	  dt <- T-mt
	  wt <- sqrt(apply(dt^2,1,sum))
	  wt <- wt/median(wt)
	  wt <- 1/((1+abs(wt/fairct))^2)
	  ngamma <- sqrt(sum(gamma^2))
	  difference <- abs(ngamma-ngammaold)/ngamma
	  wsel <- drop(wysel*wt)
          # neu BL: 17.9.09
          w0 <- which(wsel==0)
            if(length(w0) != 0)
             {
             wsel <- replace(wsel, list=w0, values=10^(-6))              
             }
          # corrected:
	  Xwsel <- Xmc[obsuse,]*sqrt(wsel)
	  ywsel <- ymc[obsuse]*sqrt(wsel)
	  #1#Xwsel <- X[obsuse,]*sqrt(wsel)
	  #1#ywsel <- y[obsuse]*sqrt(wsel)
	  loops <- loops+1
        }

########################################################################
# NEW PF:
# predictions in each segment
###b0 <- drop(t(mx)%*%b)
###ypred[-obsuse] <- drop(as.matrix(Xmc[-obsuse,]))%*%b+b0

ypred[-obsuse] <- drop(as.matrix(Xmc[-obsuse,]))%*%b
b0 <- median(y[-obsuse]-ypred[-obsuse])
ypred[-obsuse] <- ypred[-obsuse]+b0

# residuals in each segment
yresj <- y[-obsuse]-ypred[-obsuse]
SEPj[n.seg] <- sd(yresj)
#absresj=abs(yresj)
#resjlarge=quantile(absresj,1-trim)
#SEPtrimj[n.seg] <- sd(yresj[absresj<=resjlarge])
SEPtrimj[n.seg] <- sd_trim(yresj,trim=trim)
}
res <- y-ypred
########################################################################
SEPall=sd(res)
#absres=abs(res)
#reslarge=quantile(absres,1-trim)
#SEPtrim=sd(res[absres<=reslarge])
SEPtrim=sd_trim(res,trim=trim)

list(predicted=ypred,residuals=res,SEPall=SEPall,SEPtrim=SEPtrim,
	SEPj=SEPj,SEPtrimj=SEPtrimj)
} # end of prm_cvsub
#####################################################################


#####################################################################
error.bars <-
function (x, upper, lower, width = 0.02, ...)
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}
#####################################################################


#####################################################################
# start of real routine:
n=nrow(X)
p=ncol(X)

if (p>n){
        ressvd <- svd(t(X))
        X <- ressvd$v%*%diag(ressvd$d)
        n=nrow(X)
        p=ncol(X)
}

SEPall=rep(NA,a)
SEPtrim=rep(NA,a)
pred=matrix(NA,nrow=length(y),ncol=a)
SEPj=matrix(NA,nrow=segments,ncol=a)
SEPtrimj=matrix(NA,nrow=segments,ncol=a)
for (ncomp in 1:a){

print(ncomp)

  res=prm_cvsub(X=X,y=y,a=ncomp,fairct=fairct,opt=opt,subset=subset,
	segments=segments,segment.type=segment.type,trim=trim)
  SEPall[ncomp] <- res$SEPall
  SEPtrim[ncomp] <- res$SEPtrim
  pred[,ncomp] <- res$pred
  SEPj[,ncomp] <- res$SEPj
  SEPtrimj[,ncomp] <- res$SEPtrimj
}

SEPave <- apply(SEPj,2,mean)
SEPse <- apply(SEPj,2,sd)/sqrt(nrow(SEPj))
SEPtrimave <- apply(SEPtrimj,2,mean)
SEPtrimse <- apply(SEPtrimj,2,sd)/sqrt(nrow(SEPtrimj))

# optimal number of components
ind <- which.min(SEPtrimave)
fvec <- (SEPtrimave<(SEPtrimave[ind]+sdfact*SEPtrimse[ind]))
optcomp <- min((1:ind)[fvec[1:ind]])

if (plot.opt){
  # plot optimal choice of components:
  plot(1:a,SEPave,xlab="Number of components",ylab="SEP",
     ylim=range(SEPave, SEPtrimave + sdfact*SEPtrimse, SEPtrimave - sdfact*SEPtrimse),
     type="n",cex.lab=1.2)
  lines(1:a,SEPave,lty=2)
  lines(1:a,SEPtrimave,lty=1)
  points(1:a,SEPtrimave,pch=16,cex=0.5)
  error.bars(1:a, SEPtrimave + sdfact*SEPtrimse,SEPtrimave - sdfact*SEPtrimse, 
	width = 1/a,col=1)

  abline(h=SEPtrimave[ind]+sdfact*SEPtrimse[ind],lty=4)
  abline(v=optcomp,lty=4)
  legend("topright",c("SEP",paste("SEP ",trim*100,"% trimmed",sep="")),lty=c(2,1))
}


list(predicted=pred,SEPall=SEPall,SEPtrim=SEPtrim,SEPj=SEPj,SEPtrimj=SEPtrimj,optcomp=optcomp,
   SEPopt=SEPtrim[optcomp])
}

