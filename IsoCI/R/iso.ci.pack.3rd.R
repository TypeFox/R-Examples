"C.cal.logit" <- function(G.hat,C.hat){
d.logit(G.hat)*C.hat
}

"C.cal.llog" <- function(G.hat,C.hat){
d.llog(G.hat)*C.hat
}

"logit" <- function(x){
log(x/(1-x))
}

"llog" <- function(x){
log(-log(x))
}

"d.logit" <- function(x){
1/(x*(1-x))
}

"d.llog" <- function(x){
1/(x*log(x))
}


"C.cal" <- function(isofit,z0,z,d,h.opt,bw)
{
stepf <- as.stepfun(isofit)
G.hat <- stepf(z0)
g.hat <- g.density(isofit=isofit,z0=z0,z=z,d=d,h.opt=h.opt)
h.hat <- h.density(z0=z0,z=z,bw=bw)

C.hat <- ((4*g.hat*G.hat*(1-G.hat))/h.hat)^(1/3)
C.hat
}

"g.density" <- function(isofit,z0,z,d,h.opt)
{
n <- length(d) # the number of obs.
G.hat <- isofit$yf
stepf <- as.stepfun(isofit)

G.z0.hat <- stepf(z0)
G.diff <- c(G.hat[1],diff(G.hat))
g <- dnorm((z0-z)/h.opt)/h.opt * G.diff
g <- sum(g)
g
}


"h.density" <- function(z0,z,bw)
{
#require(KernSmooth)
#h <- dpik(z)
n <- length(z)
Density <- dnorm((z0-z)/bw)/bw
Density <- mean(Density)
Density
}

"iso.ci.transform.unit" <- function(isofit,z0,z,d,alpha,h.opt,bw){
alpha.set <- c(0.2,0.1,0.05,0.02,0.01,0.002)
if(!(alpha %in% alpha.set)) stop("invalid alpha")

if(alpha==0.2) q.alpha <- 0.664235
if(alpha==0.1) q.alpha <- 0.845081
if(alpha==0.05) q.alpha <- 0.998181
if(alpha==0.02) q.alpha <- 1.171530
if(alpha==0.01) q.alpha <- 1.286659
if(alpha==0.002) q.alpha <- 1.516664

n <- length(z)
stepf <- as.stepfun(isofit)
G.hat <- stepf(z0)
C.hat <- C.cal(isofit=isofit,z0=z0,z=z,d=d,h.opt=h.opt,bw=bw)
C.hat.logit <- C.cal.logit(G.hat=G.hat,C.hat=C.hat)
C.hat.llog <- C.cal.logit(G.hat=G.hat,C.hat=C.hat)
lhm.un <- G.hat - (n^(-1/3)) * C.hat * q.alpha
rhm.un <- G.hat + (n^(-1/3)) * C.hat * q.alpha
lhm.1 <- logit(G.hat) - (n^(-1/3)) * C.hat.logit * q.alpha
rhm.1 <- logit(G.hat) + (n^(-1/3)) * C.hat.logit * q.alpha
lhm.2 <- llog(G.hat) - (n^(-1/3)) * C.hat.llog * q.alpha
rhm.2 <- llog(G.hat) + (n^(-1/3)) * C.hat.llog * q.alpha
lhm.logit <- exp(lhm.1)/(1+exp(lhm.1)) 
rhm.logit <- exp(rhm.1)/(1+exp(rhm.1))
lhm.llog <- exp(-exp(rhm.2))
rhm.llog <- exp(-exp(lhm.2))
if(G.hat == 0){lhm.un <- rhm.un <- lhm.logit <- rhm.logit <- lhm.llog <- rhm.llog <- 0}
if(G.hat == 1){lhm.un <- rhm.un <- lhm.logit <- rhm.logit <- lhm.llog <- rhm.llog <- 1}
cbind(G.hat,lhm.un,rhm.un,lhm.logit,rhm.logit,lhm.llog,rhm.llog)
}

"iso.bt.unit.sub" <- function(dd.vec,z0){
n2 <- length(dd.vec)
n <- n2/2
z <- dd.vec[1:n]
d <- dd.vec[-c(1:n)]
isofit.boots <- isoreg(z,d)
stepfn <- as.stepfun(isofit.boots)
G.hat <- stepfn(z0)
}

"iso.bt.unit" <- function(dd,z0,alpha){
G.hat.vec <- apply(dd,2,iso.bt.unit.sub,z0=z0)
q <- quantile(G.hat.vec,c(alpha/2,1-alpha/2))
q
}

"iso.bt.wald.unit" <- function(dd,z0){
G.hat.vec <- apply(dd,2,iso.bt.unit.sub,z0=z0)
var(G.hat.vec)
}



"bandwidth.choose.unit" <- function(h.val,z,d){
n <- length(z)
id <- 1:n
g.hat.cv <- mapply(bandwidth.choose.cv,id=id,MoreArgs=list(z=z,d=d,h.val=h.val))
sum(log(g.hat.cv))
}

"bandwidth.choose.cv" <- function(id,z,d,h.val){
z.test <- z[id]
z.train <- z[-id]
d.train <- d[-id]
isofit <- isoreg(z.train,d.train)
g.hat <- g.density(isofit=isofit,z0=z.test,z=z.train,d=d.train,h.opt=h.val)
g.hat[g.hat == 0] <- 1e-20
g.hat
}


bandwidth.choose <- function(h.set,z,d){
LCV <- mapply(bandwidth.choose.unit,h.val=h.set,MoreArgs=list(z=z,d=d))
result.table <- cbind(h.set,LCV)
h.opt <- min(h.set[LCV==max(LCV)])
list(h.opt=h.opt,result.table=result.table)
}


iso.ci <- function(z,d,alpha=0.05,h.opt=0.3,
	nboots=500,method="wald.tr",seed=1253)
{

if(any(is.na(z)) || any(is.na(d))) 
	stop("missing values not allowed")

if(length(z) != length(d))
	stop("lengths of z and d must be the same.")

if(!(method %in% c("wald.tr","bt","bt.wald"))) stop("invalid method")

n <- length(z)
o <-  order(z)
z <- z[o]
d <- d[o]

isofit <- isoreg(z,d)
yf <- isofit$yf



if(method == "wald.tr")
{

alpha.set <- c(0.2,0.1,0.05,0.02,0.01,0.002)
if(!(alpha %in% alpha.set)) stop("invalid alpha for transformed CI")

if(alpha==0.2) q.alpha <- 0.664235
if(alpha==0.1) q.alpha <- 0.845081
if(alpha==0.05) q.alpha <- 0.998181
if(alpha==0.02) q.alpha <- 1.171530
if(alpha==0.01) q.alpha <- 1.286659
if(alpha==0.002) q.alpha <- 1.516664

#require(KernSmooth)
bw <- dpik(z)

output <- mapply(iso.ci.transform.unit,z0=z,
		MoreArgs=list(isofit=isofit,z=z,d=d,alpha=alpha,h.opt=h.opt,bw=bw))
output1 <- t(output)
output2 <- cbind(z,output1)
lhm <- output2[,c(3,5,7)]
rhm <- output2[,c(4,6,8)]
lhm.bd <- ifelse(lhm <0, 0, ifelse(lhm > 1, 1, lhm)) 
rhm.bd <- ifelse(rhm <0, 0, ifelse(rhm > 1, 1, rhm))
output2[,c(3,5,7)] <- lhm.bd
output2[,c(4,6,8)] <- rhm.bd
colnames(output2) <- c("z","yf","wald.lhm","wald.rhm",
	"logit.lhm","logit.rhm","llog.lhm","llog.rhm")
return(output2)
} else if(method == "bt"){
set.seed(seed)
dd <- matrix(0,2*n,nboots)
	for(i in 1:nboots){
	ss <- sample(1:n,n,TRUE)
	z.ss <- z[ss]
	d.ss <- d[ss]
	o.ss <- order(z.ss)
	z.ss <- z.ss[o.ss]
	d.ss <- d.ss[o.ss]
	dd[,i] <- c(z.ss,d.ss)
	}
output <- mapply(iso.bt.unit,z0=z,MoreArgs=list(dd=dd,alpha=alpha))
output1 <- t(output)
output2 <- cbind(z,yf,output1)
colnames(output2) <- c("z","yf","nbt.lhm","nbt.rhm")
return(output2) 
} else if(method == "bt.wald"){
alpha.set <- c(0.2,0.1,0.05,0.02,0.01,0.002)
if(!(alpha %in% alpha.set)) stop("invalid alpha for transformed CI")

if(alpha==0.2) q.alpha <- 0.664235
if(alpha==0.1) q.alpha <- 0.845081
if(alpha==0.05) q.alpha <- 0.998181
if(alpha==0.02) q.alpha <- 1.171530
if(alpha==0.01) q.alpha <- 1.286659
if(alpha==0.002) q.alpha <- 1.516664

set.seed(seed)
dd <- matrix(0,2*n,nboots)
	for(i in 1:nboots){
	ss <- sample(1:n,n,TRUE)
	z.ss <- z[ss]
	d.ss <- d[ss]
	o.ss <- order(z.ss)
	z.ss <- z.ss[o.ss]
	d.ss <- d.ss[o.ss]
	dd[,i] <- c(z.ss,d.ss)
	}
var.Ghat <- mapply(iso.bt.wald.unit,z0=z,MoreArgs=list(dd=dd))
C.hat <- (n^(1/3))*sqrt(var.Ghat/0.2636)
G.hat <- yf
C.hat.logit <- C.cal.logit(G.hat=G.hat,C.hat=C.hat)
C.hat.llog <- C.cal.logit(G.hat=G.hat,C.hat=C.hat)
lhm.un <- G.hat - (n^(-1/3)) * C.hat * q.alpha
rhm.un <- G.hat + (n^(-1/3)) * C.hat * q.alpha
lhm.1 <- logit(G.hat) - (n^(-1/3)) * C.hat.logit * q.alpha
rhm.1 <- logit(G.hat) + (n^(-1/3)) * C.hat.logit * q.alpha
lhm.2 <- llog(G.hat) - (n^(-1/3)) * C.hat.llog * q.alpha
rhm.2 <- llog(G.hat) + (n^(-1/3)) * C.hat.llog * q.alpha
lhm.logit <- exp(lhm.1)/(1+exp(lhm.1)) 
rhm.logit <- exp(rhm.1)/(1+exp(rhm.1))
lhm.llog <- exp(-exp(rhm.2))
rhm.llog <- exp(-exp(lhm.2))
output2 <- cbind(z,G.hat,lhm.un,rhm.un,lhm.logit,rhm.logit,lhm.llog,rhm.llog) 

for(i in 1:n){
	if(G.hat[i]==0) {output2[i,3:8] <- 0}
	else if(G.hat[i]==1) {output2[i,3:8] <- 1}
	else {output2[i,3:8] <- output2[i,3:8]}
}

lhm <- output2[,c(3,5,7)]
rhm <- output2[,c(4,6,8)]
lhm.bd <- ifelse(lhm <0, 0, ifelse(lhm > 1, 1, lhm)) 
rhm.bd <- ifelse(rhm <0, 0, ifelse(rhm > 1, 1, rhm))
output2[,c(3,5,7)] <- lhm.bd
output2[,c(4,6,8)] <- rhm.bd

colnames(output2) <- c("z","yf","bt.wald.lhm","bt.wald.rhm",
	"bt.logit.lhm","bt.logit.rhm","bt.llog.lhm","bt.llog.rhm")
output2
}
}

