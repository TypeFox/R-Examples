`tlrt` <-
function (y, p, d=1, transform = "no",  a = 0.25, b = 0.75,...)
{
makedata <-
function (dataf, p1, p2, d, is.constant1 = TRUE, is.constant2 = TRUE, 
    thd.by.phase = FALSE) 
{
    n <- dim(dataf)[1]
    nseries <- dim(dataf)[2]
    start <- max(c(p1, p2, d)) + 1
    p <- max(c(p1, p2))
    xy <- NULL
    for (i in (1:nseries)) xy <- rbind(xy, cbind(makexy(dataf[, 
        i], p, start, d, thd.by.phase = thd.by.phase), i))
    xy <- dna(xy)
    xy <- xy[s <- order(xy[, p + 2]), ]
    xy1 <- setxy(xy, p, p1, nseries, is.constant1)
    xy2 <- setxy(xy, p, p2, nseries, is.constant2)
    list(xy1 = xy1, xy2 = xy2, sort.list = s)
}

makexy <-
function (x, p, start, d, thd.by.phase = FALSE) 
{
    n <- length(x)
    xy <- NULL
    for (i in (1:p)) xy <- cbind(xy, x[(start - i):(n - i)])
    if (thd.by.phase) 
        xy <- cbind(xy, x[start:n], x[(start - 1):(n - 1)] - 
            x[(start - 2):(n - 2)])
    else xy <- cbind(xy, x[start:n], x[(start - d):(n - d)])
    xy
}

setxy<-
function (old.xy, p, p1, nseries, is.coefficient = TRUE) 
{
    if (p1 >= 1) 
        s <- 1:p1
    else s <- NULL
    n <- dim(old.xy)[1]
    new.xy <- old.xy[, c(s, (p + 1):(p + 3)), drop = FALSE]
    temp <- new.xy[, s, drop = FALSE]
    if (is.coefficient) 
        new.xy <- cbind(1, new.xy)
    if (is.coefficient) 
        temp <- cbind(rep(1, n), temp)
    if (nseries == 1) 
        return(new.xy)
    for (i in rev(2:nseries)) {
        select <- old.xy[, p + 3] == i
        zero <- 0 * temp
        zero[select, ] <- temp[select, ]
        new.xy <- cbind(zero, new.xy)
    }
    new.xy
}

findstart <-
function (x, nseries, indexid, p) 
{
    m <- dim(x)[1]
    amax <- 0
    for (i in (1:nseries)) {
        amax <- max(amax, (1:m)[(cumsum(x[, indexid] == i) == 
            p)])
    }
    amax
}

dna <-
function (x) 
{
    x[!apply(x, 1, any.na), ]
}

any.na <-
function (x) 
{
    any(x == "NA")
}

revm <-
function (m) 
{
    apply(m, 2, rev)
}
p.value.tlrt <-
function(y,a=0.25,b=0.75,p=0){
# This function computes the approximate p-value of the threshold likelihood
# ratio test as reported in Chan (1990).
# input:
# y=the likelihood ratio test statistic
# a,b = the test assumes that the threshold is searched from the a*100 to
#       b*100 percentiles of the data
# p= autoregressive order 
# The delay d is assumed to be <=p
#
# output:
# p-value of the likelihood ratio test for threshold linearity
#
#
t1=function(x){F=pnorm(x)
0.5*log(F/(1-F))
}
lower=qnorm(a)
upper=qnorm(b)
#if(p==0) temp=integrate(t1,lower,upper)$value
if(p==0) temp=t1(upper)-t1(lower)
if (p>0) {
tp1=function(x){
F=pnorm(x)
f=dnorm(x)
b=2*F-x*f
c=F*(F-x*f)-f*f
root=0.5*(b+sqrt(b*b-4*c))
0.5*log(root/(1-root))
}
tp2=function(x){
F=pnorm(x)
f=dnorm(x)
b=2*F-x*f
c=F*(F-x*f)-f*f
root=0.5*(b-sqrt(b*b-4*c))
0.5*log(root/(1-root))
}

#temp=(p-1)*integrate(t1,lower,upper)$value+integrate(tp1,lower,upper)$value+
#integrate(tp2,lower,upper)$value
temp=(p-1)*(t1(upper)-t1(lower))+tp1(upper)-tp1(lower)+tp2(upper)-tp2(lower)
}
if(p>0) return(1-exp(-2*dchisq(y,df=p+1)*(y/(p+1)-1)*temp)) else {
z=sqrt(y)
alpha=temp
return(sqrt(2/pi)*exp(-y/2)*(temp*(z-1/z)+1/z))
}
}

# Conduct the likelihood ratio test for threshold nonlinearity
#
dataf=y
if(missing(p)) p=ar(y,...)$order 

p1=p
p2=p
                if (!is.matrix(dataf)) {
        temp <- cbind(dataf, dataf)
        dataf <- temp[, 1, drop = FALSE]
    }
    dataf <- switch(transform, log = log(dataf), log10 = log10(dataf), 
        sqrt = sqrt(dataf), no = dataf)
          res <- makedata(dataf, p1, p2, d)
    xy1 <- res$xy1
    xy2 <- res$xy2
    sort.l <- res$sort.list
    m <- dim(xy1)[1]
    q1 <- dim(xy1)[2]
    rss1=rep(10^10,m)
    rss2=rep(10^10,m)
    s <- (q1 - 3 - p1):(q1 - 3)
    temp <- xy1[, s]
    xy1 <- cbind(temp, xy1[, -s])
             lbound <- sum(dataf[dataf != "NA"] == min(dataf[dataf != 
        "NA"]))
    ubound <- sum(dataf[dataf != "NA"] == max(dataf[dataf != 
        "NA"]))
    i1 <- max(c(q1 - 3, lbound + 1,  p + 1, d, findstart(xy1, 
        1, q1, p1 + 2)))
    i1 <- max(i1, floor(a * m))
    i2 <- m - max(c(q1 - 3, ubound + 1, p + 1, d, findstart(revm(xy1), 
        1, q1, p1 + 2))) - 1
    i2 <- min(i2, ceiling(b * m))
#    i1=ceiling(a*m)
#    i2=floor(b*m) 
truea=i1/m
trueb=i2/m
    s <- -((q1 - 1):q1)
    R <- qr.R(qr(xy1[1:i1, s]))
        posy <- q1 - 2
        rss1[i1] <- (R[posy, posy])^2
                for (i in ((i1 + 1):i2)) {
            R <- qr.R(qr(rbind(R, xy1[i, s])))
            rss1[i] <- (R[posy, posy])^2
            }
        s <- -((q1 - 1):q1)
        posy <- q1 - 2
        R <- qr.R(qr(xy1[(i2 + 1):m, s]))
        rss2[i2] <- (R[posy, posy])^2
                for (i in rev(i1:(i2 - 1))) {
            R <- qr.R(qr(rbind(R, xy1[i + 1, s])))
            rss2[i] <- (R[posy, posy])^2
        }
        rss <- rss1 + rss2
        rss.H1=min(rss)
        R=qr.R(qr(xy1[,s]))
        rss.H0=(R[posy,posy])^2
        test.stat=m*(rss.H0-rss.H1)/rss.H1
list(percentiles=signif(c(truea,trueb)*100,3),test.statistic=test.stat,
p.value=p.value.tlrt(test.stat,a=truea,b=trueb,p=p))
}

