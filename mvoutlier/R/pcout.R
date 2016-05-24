"pcout" <-
function(x,makeplot=FALSE,
      explvar=0.99,crit.M1=1/3,crit.c1=2.5,crit.M2=1/4,crit.c2=0.99,
      cs=0.25,outbound=0.25, ...){

#################################################################
p=ncol(x)
n=nrow(x)

x.mad=apply(x,2,mad)
if (any(x.mad==0))
  stop("More than 50% equal values in one or more variables!")

#################################################################
# PHASE 1:

# Step 1: robustly sphere the data:
x.sc <- scale(x,apply(x,2,median),x.mad)

# Step 2: PC decomposition; compute p*, robustly sphere:
##x.wcov <- cov(x.sc)
##x.eig <- eigen(x.wcov)
##p <- ncol(x)
##p1 <- (1:p)[(cumsum(x.eig$val)/sum(x.eig$val)>explvar)][1]

# or faster:
x.svd <- svd(scale(x.sc,TRUE,FALSE))
a <- x.svd$d^2/(n-1)
p1 <- (1:p)[(cumsum(a)/sum(a)>explvar)][1]

x.pc <- x.sc%*%x.svd$v[,1:p1]
xpc.sc <- scale(x.pc,apply(x.pc,2,median),apply(x.pc,2,mad))

# Step 3: compute robust kurtosis weights, transform to distances:
wp <- abs(apply(xpc.sc^4,2,mean)-3)

xpcw.sc <- xpc.sc%*%diag(wp/sum(wp))
xpc.norm <- sqrt(apply(xpcw.sc^2,1,sum))
x.dist1 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)

# Step 4: determine weights according to translated biweight:
M1 <- quantile(x.dist1,crit.M1)
const1 <- median(x.dist1)+crit.c1*mad(x.dist1)
w1 <- (1-((x.dist1-M1)/(const1-M1))^2)^2
w1[x.dist1<M1] <- 1
w1[x.dist1>const1] <- 0

#################################################################
# PHASE 2:

# Step 5: compute Euclidean norms of PCs and their distances:
xpc.norm <- sqrt(apply(xpc.sc^2,1,sum))
x.dist2 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)

# Step 6: determine weight according to translated biweight:
M2 <- sqrt(qchisq(crit.M2,p1))
const2 <- sqrt(qchisq(crit.c2,p1))
w2 <- (1-((x.dist2-M2)/(const2-M2))^2)^2
w2[x.dist2<M2] <- 1
w2[x.dist2>const2] <- 0

#################################################################
# Combine PHASE1 and PHASE 2: compute final weights:
wfinal <- (w1+cs)*(w2+cs)/((1+cs)^2)
wfinal01 <- round(wfinal+0.5-outbound)

#################################################################
# Generate plot:
if (makeplot){
  op <- par(mfrow=c(3,2), mar=c(4,4,2,2))
  on.exit(par(op))

  # location outliers
  plot(x.dist1,xlab="Index",ylab="Distance (location)", ...)
  abline(h=const1)
  abline(h=M1,lty=2)
  plot(w1,xlab="Index",ylab="Weight (location)",ylim=c(0,1), ...)
  abline(h=0)
  abline(h=1,lty=2)

  # scatter outliers
  plot(x.dist2,xlab="Index",ylab="Distance (scatter)", ...)
  abline(h=const2)
  abline(h=M2,lty=2)
  plot(w2,xlab="Index",ylab="Weight (scatter)",ylim=c(0,1), ...)
  abline(h=0)
  abline(h=1,lty=2)

  # combined weights
  plot(wfinal,xlab="Index",ylab="Weight (combined)",ylim=c(0,1), ...)
  abline(h=cs)
  plot(wfinal01,xlab="Index",ylab="Final 0/1 weight",ylim=c(0,1), ...)
}

list(wfinal01=wfinal01,wfinal=wfinal,wloc=w1,wscat=w2,
     x.dist1=x.dist1,x.dist2=x.dist2,M1=M1,const1=const1,M2=M2,const2=const2)
}

