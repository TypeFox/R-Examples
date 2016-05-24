pamr.score.to.class1 <- function (x, scores, cutoff=2, n.class=2) 
{
        x.sml <- x[abs(scores)>cutoff,]
        out <- kmeans2(t(x.sml), centers=n.class)$cluster 
        return(out)
}
pamr.score.to.class2 <- function (x, scores, cutoff=2, n.pc=1, n.class=2) 
{
        x.sml <- x[abs(scores)>cutoff,]
        v <- svd(x.sml)$v
        v.sml <- v[,1:n.pc]
        out <- kmeans2(v.sml, centers=n.class)$cluster 
        return(out)
}

pamr.surv.to.class2 <- function (y, icens, cutoffs=NULL, n.class=NULL,  class.names=NULL, newy=y, newic=icens)

# Splits patients into classes based on their survival times
# The user can either specify the number of classes or the survival
# time cutoffs.
#
# y - vector of survival times
# icens - censoring indicator
# cutoffs - survival time cutoffs
# n.class - number of classes to create
# class.names - optional vector of names for each class
{
 #       require(survival)
        if (is.null(cutoffs) & is.null(n.class)) {
                stop("Must specify either cutoffs or n.class")
        }
        if (!is.null(cutoffs) & !is.null(n.class)) {
                stop("Can't have both cutoffs and n.class specified")
        }
        data.sfit <- survfit(Surv(y,icens)~1)
        if (!is.null(cutoffs)) {
                if (is.null(class.names)) {
                        class.names <- 1:(length(cutoffs)+1)
                }
                cur.mat <- gen.y.mat2(list(y=y, icens=icens), cutoffs, class.names,                                              newdata=list(y=newy, icens=newic))
        }
        else {
                if (n.class==1) {
                        stop("Must have at least two classes")
                }
                if (is.null(class.names)) {
                        class.names <- 1:n.class
                }
                cur.quantiles <- seq(from=0, to=1, length=n.class+1)
                cur.quantiles <- cur.quantiles[2:n.class]
                cutoffs <- quantile(y[icens==1], cur.quantiles)
                cur.mat <- gen.y.mat2(list(y=y, icens=icens), cutoffs, class.names,
                                newdata=list(y=newy, icens=newic))
        }
        mle.classes <- apply(cur.mat, 1, get.mle.class)
         return(list(class=as.numeric(mle.classes), prob=cur.mat, cutoffs=cutoffs))
}
gen.y.mat2 <- function(surv.data, cutoffs, class.names=NULL, newdata=surv.data)
# Calculates the probability that a given patient belongs to a given
# class.  Returns a matrix where entry (i,j) is the probability that
# patient i belongs to class j.  The function calculates the
# probability that a given patient dies between two given cutoffs,
# and uses this information to calculate the probability that
# a patient with a censored survival time died in a given interval.
{
         data.sfit <- survfit(Surv(surv.data$y,surv.data$icens)~1)
         surv.ndx <- find.surv.ndx(cutoffs, data.sfit$time)
         surv.probs <- c(0, 1-data.sfit$surv[surv.ndx],1)
         surv.probs <- c(rep(0, sum((surv.ndx==0))), surv.probs)
         cutoffs <- c((min(surv.data$y)-1), cutoffs, (max(surv.data$y)+1))
         y.size <- length(cutoffs)
         y.mat <- matrix(0,nrow=length(newdata$y), ncol=(y.size-1))
         for (i in 2:y.size) {
                 cur.int.prob <- surv.probs[i] - surv.probs[i-1]
                 y.mat[((newdata$y<=cutoffs[i])&(newdata$y>cutoffs[i-1])&
                         (newdata$icens==1)),i-1] <- 1
                 which.x <- ((newdata$icens==0)&(newdata$y<=cutoffs[i-1]))
                 if (sum(which.x)>0) {
                         which.x.vals <- newdata$y[which.x]
                         surv.ndx <- find.surv.ndx(which.x.vals,
                                 data.sfit$time)
                         y.mat[which.x,i-1][surv.ndx==0] <- cur.int.prob
                         y.mat[which.x,i-1][surv.ndx!=0] <- cur.int.prob /
                                 data.sfit$surv[surv.ndx]
                 }
                 which.x <- ((newdata$icens==0)&(newdata$y>cutoffs[i-1])&
                         (newdata$y<=cutoffs[i]))
                 if (sum(which.x>0)) {
                         which.x.vals <- newdata$y[which.x]
                         surv.ndx <- find.surv.ndx(which.x.vals,
                                 data.sfit$time)
                         y.mat[which.x,i-1][surv.ndx==0] <- surv.probs[i]
                         y.mat[which.x,i-1][surv.ndx!=0] <- 1 -
                                 (1 - surv.probs[i]) / data.sfit$surv[surv.ndx]
                 }
         }
         if (!is.null(class.names)) {
                 y.mat <- as.data.frame(y.mat)
                 names(y.mat) <- class.names
                 y.mat <- as.matrix(y.mat)
         }
         y.mat
}
  


get.surv.q <- function(surv.obj, quantile) 
{
    ndx <- sum(surv.obj$surv > quantile)
    if (ndx==0)
        return(0)
    else
        return(surv.obj$time[ndx])
}
find.surv.ndx <- function(newtimes, oldtimes) 
{
	first <- apply(as.matrix(newtimes), 1, function(e1,e2) (e1>=e2), e2=oldtimes)
	as.vector(apply(first, 2, sum))
}
get.mle.class <- function(y.row) 
{
	i <- 1+sum((max(y.row)>cummax(y.row)))
	if (!is.null(names(y.row)[i])) {
		return(names(y.row)[i])
	}
	else return(i)
}

kmeans2 <- function(x, ..., n.rep=10) 
# Performs k-means clustering multiple times from different starting
# points
{
#        require(cluster)
        wss <- Inf
        for (i in 1:n.rep) {
                cur.fit <- kmeans(x, ...)
                if (sum(cur.fit$withinss) < wss) {
                        fit <- cur.fit
                        wss <- sum(fit$withinss)
                }
        }
        return(fit)
}
sam.func <- function (x, y, fudge=median(sd)) 
{
        y.l <- levels(y)
        x.n1 <- sum(y==y.l[1])
        x.n2 <- sum(y==y.l[2])
        x.1 <- x[,(y==y.l[1])]
        x.2 <- x[,(y==y.l[2])]
        x.bar1 <- apply(x.1, 1, mean)
        x.bar2 <- apply(x.2, 1, mean)
        x.var1 <- apply(x.1, 1, var) * (x.n1-1)
        x.var2 <- apply(x.2, 1, var) * (x.n2-1)
        sd <- sqrt(((1/x.n1 + 1/x.n2)/(x.n1+x.n2-2)) * (x.var1+x.var2))
        numer <- (x.bar1 - x.bar2)
        tt <- numer / (sd+fudge)
        return(list(tt=tt, numer=numer, sd=sd))
}
"cox.func"<-
function(x, y, icens, fudge = 0)
{
	scor <- coxscor(x, y, icens)$scor
	sd <- sqrt(coxvar(x, y, icens))
	tt <- scor/(sd + fudge)
	return(tt, numer = scor, sd)
}
"coxscor"<-
function(x, y, ic, offset = rep(0, length(y)))
{
# computes cox scor function for rows of nx by n matrix  x
	n <- length(y)
	nx <- nrow(x)
	yy <- y + (ic == 0) * (1.0000000000000001e-05)
	otag <- order(yy)
	y <- y[otag]
	ic <- ic[otag]
	x <- x[, otag, drop = F]
	offset <- offset[otag]	
	#compute  unique failure times, d=# of deaths at each failure time, 
#dd= expanded version of d to length n, s=sum of covariates at each
# failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time
	a <- coxstuff(x, y, ic, offset = offset)
	nf <- a$nf
	fail.times <- a$fail.times
	s <- a$s
	d <- a$d
	dd <- a$dd
	nn <- a$nn
	nno <- a$nno
	w <- rep(0, nx)
	for(i in (1:nf)) {
		w <- w + s[, i]
		for(j in (1:n)[y >= fail.times[i]]) {
			w <- w - (d[i] * x[, j, drop = F] * safe.exp(offset[j]))/nno[
				i]
		}
	}
	return(scor = w, coxstuff.obs = a)
}
"coxstuff"<-
function(x, y, ic, offset = rep(0, length(y)))
{
	fail.times <- unique(y[ic == 1])
	nf <- length(fail.times)
	n <- length(y)
	nn <- rep(0, nf)
	nno <- rep(0, nf)
	for(i in 1:nf) {
		nn[i] <- sum(y >= fail.times[i])
		nno[i] <- sum(safe.exp(offset)[y >= fail.times[i]])
	}
	s <- matrix(0, ncol = nf, nrow = nrow(x))
	d <- rep(0, nf)
	for(i in 1:nf) {
		o <- (1:n)[(y == fail.times[i]) & (ic == 1)]
		d[i] <- length(o)
		s[, i] <- apply(x[, o, drop = F], 1, sum)
	}
#expand d out to a vector of length n
	dd <- rep(0, n)
	for(j in 1:nf) {
		dd[(y == fail.times[j]) & (ic == 1)] <- d[j]
	}
	return(fail.times, s, d, dd, nf, nn, nno)
}
"coxvar"<-
function(x, y, ic, offset = rep(0, length(y)), coxstuff.obj = NULL)
{
# computes information elements (var) for cox
# x is nx by n matrix of expression  values
	nx <- nrow(x)
	n <- length(y)
	yy <- y + (ic == 0) * (9.9999999999999995e-07)
	otag <- order(yy)
	y <- y[otag]
	ic <- ic[otag]
	x <- x[, otag, drop = F]
	offset <- offset[otag]
	if(is.null(coxstuff.obj)) {
		coxstuff.obj <- coxstuff(x, y, ic, offset = offset)
	}
	nf <- coxstuff.obj$nf
	fail.times <- coxstuff.obj$fail.times
	s <- coxstuff.obj$s
	d <- coxstuff.obj$d
	dd <- coxstuff.obj$dd
	nn <- coxstuff.obj$nn
	nno <- coxstuff.obj$nno
	w <- rep(0, nx)
	for(i in 1:nf) {
		sx <- rep(0, nx)
		s <- rep(0, nx)
		ii <- (1:n)[y >= fail.times[i]]
		for(j in ii) {
			sx <- sx + (x[, j] * safe.exp(offset[j]))/nno[i]
		}
		for(j in ii) {
			s <- s + (x[, j]^2 * safe.exp(offset[j]))/nno[i]
		}
		w <- w + d[i] * (s - sx * sx)
	}
	return(w)
}
pamr.confusion.survival <-
  
 function(fit, survival.time,censoring.status, yhat){
   
# computes confusion matrix for (survival.time,censoring) outcome
# based on fit object "fit" and class predictions "yhat"
# soft response probabilities for (survival.time,censoring) are first estimated
#  using Kaplan-Meier method applied to training data
   
  n.class<- fit$ngroup.survival
  
true <- pamr.surv.to.class2(fit$survival.time, fit$censoring.status,newy=survival.time, newic=censoring.status, n.class=n.class)$prob

# use of temp below is to ensure that a categeory is not dropped
temp <- c(yhat, 1:n.class)
Yhat <- model.matrix( ~ factor(temp) - 1,data = list(y = temp))
Yhat <- Yhat[1:length(yhat),]


ytemp<- apply(true,1,which.is.max)
 temp <- c(yhat,names(table(ytemp)))
   nams <- names(table(temp))

      tt <- matrix(NA,nrow=length(fit$prior),ncol=length(fit$prior))
  
         for(i in 1:length(fit$prior)){
           for(j in 1:length(fit$prior)){
                 tt[i,j] <- sum(true[,i]*Yhat[,j])
               }}
     dimnames(tt) <- list(names(table(ytemp)),nams)

   tt1 <- tt
    diag(tt1) <- 0
    tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
    dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"

return(tt)
}
pamr.plotsurvival <- function(group, survival.time, censoring.status){
  # plots Kaplan-Meier curves stratified by "group"
#  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status)~as.factor(group))
  junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))      
  pv <- 1-pchisq(2*(junk2$loglik[2]-junk2$loglik[1]),df=n.class-1)

 plot(junk, col=2:(2+n.class-1) ,xlab= "Time", ylab="Probability of survival")
  
 legend(.01*max(survival.time),.2, col=2:(2+n.class-1), lty=rep(1,n.class),
        legend=as.character(1:n.class))
 text(0.01 * max(survival.time), .25, paste("pvalue=",as.character(round(pv,4))))
return()
}
 pamr.pvalue.survival <- function(group, survival.time, censoring.status,
                ngroup.survival)
        {
                temp <- coxph(Surv(survival.time, censoring.status) ~ as.factor(
                        group))
                loglik <- 2 * (temp$loglik[2] - temp$loglik[1])
                return(1 - pchisq(loglik, ngroup.survival - 1))
        }
order.class.survival <- function(a, survival.time, censoring.status){
#require(survival)

#orders the classes specified in "a" by median survival time, from
#shortest to longest

med <- rep(NA,length(table(a)))
for(i in 1:length(table(a))){
    o <- a==i
    aa <- survfit(Surv(survival.time[o],censoring.status[o]))
    med[i] <- approx(aa$surv,aa$time, xout=.5, method="constant")$y
    }
aa <- rep(NA,length(a))
for(i in 1:length(table(a))){
  aa[a==i] <- rank(med)[i]
}
return(aa)
}

pamr.test.errors.surv.compute <- function(proby, yhat) {
## computes confusion matrix, class-wise error rate and overall error rate
## rows of confusion matrix refer to true classes; columns to predicted classes
##
## yhat is class prediction, proby is matrix of "true" soft class probabilities.
## 
  true <- proby
  ytemp <- apply(true, 1, which.is.max)
  temp <- c(yhat, names(table(ytemp)))
  nams <- names(table(temp))
  Yhat <- model.matrix(~factor(temp) - 1, data = list(y = temp))
  Yhat <- Yhat[1:length(yhat), ]
  nc <- ncol(proby)
  tt <- matrix(NA, nrow = nc, ncol = nc)
  for (i in 1:nc) {
    for (j in 1:nc) {
      tt[i, j] <- sum(true[, i] * Yhat[, j])
    }
  }
  dimnames(tt) <- list(names(table(ytemp)), nams)
  
  tt1 <- tt
  diag(tt1) <- 0
  tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
  dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
  error <- sum(tt1)/sum(tt)
  return(list(confusion=tt,error=error))
}

safe.exp=function(x){
 xx=sign(x)*pmin(abs(x),500)
 return(exp(xx))
}

