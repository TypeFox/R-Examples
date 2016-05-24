boot.ratio.sd.bca <-
function(x, y, stacked = TRUE, variable = NULL, null.hyp = NULL,
                              alternative = c("two.sided","less","greater"), conf.level = 0.95,
                              type = NULL, R = 9999)
{
# require(boot)
if (stacked)
  {
   variable <- all.names(substitute(x))
   if (length(variable)>1) variable <- variable[[3]]
   y <- as.factor(y)
   l <- levels(y)
   pop.1 <- l[1]
   pop.2 <- l[2]
   u <- x[y==l[1]]
   v <- x[y==l[2]]
   x <- u
   y <- v
  } 
else
  {
   pop.1 <- all.names(substitute(x))
   if (length(pop.1)>1) pop.1 <- pop.1[[3]]
   pop.2 <- all.names(substitute(y))
   if (length(pop.2)>1) pop.2 <- pop.2[[3]]
  }
test <- !is.null(null.hyp)
m <- length(x)
n <- length(y)
obs.1 <- sd(x)
obs.2 <- sd(y)
obs <- obs.1/obs.2
sd.fun <- function(w,i) sd(w[i])
boot.sd.x <- boot(x,sd.fun,R=R)$t
boot.sd.y <- boot(y,sd.fun,R=R)$t
z <- c(boot.sd.x/boot.sd.y)
bias <- signif(mean(z)-obs,digits=3)
percent.bias <- signif(100*abs(bias/obs),digits=3)
cl <- paste(100*conf.level,"%",sep="")
if (identical(alternative,c("two.sided","less","greater"))) alternative <- "two.sided"
z1 <- mean(z < obs)
z0 <- qnorm(z1)
u <- sd(x)
v <- sd(y)
d <- vector(length=m+n)
for (i in 1:m) {d[i] <- sd(x[-i])-v}
for (j in 1:n) {d[m+j] <- u-sd(y[-j])}
mean.d <- mean(d)
a0 <- -sum((d-mean.d)^3)/(6*(sum((d-mean.d)^2))^(3/2))
# FOR THE HYPOTHESIS TEST
if(test)
{
 rtp <- (sum((z < null.hyp))+(sum((z == null.hyp))+1)/2)/(R+1)
 b0 <- qnorm(rtp)
 c0 <- ((2+a0*b0-a0*z0)*z0-b0)/(1+a0*b0-a0*z0)
 p0 <- pnorm(c0) # P-value for a left-tailed test
 tc <- c("two.sided","less","greater")
 pc <- c(2*min(p0,1-p0),p0,1-p0)
 p <- signif(pc[tc==alternative],digits=3)
 pv <- c((p>=0.001)&(p<=0.999),(p<0.001),(p>0.999))
 pt <- c(p,"P < 0.001","P > 0.999")
 p.value <- pt[pv]
 ac <- c("not-equal","less-than","greater-than")
 alt <- ac[tc==alternative]
}
# FOR THE CONFIDENCE INTERVAL
tci <- c("two.sided","greater","less")
ti <- c("two-sided","lower-bound","upper-bound")
if (is.null(type)) type <- ti[tci==alternative]
a <- if (identical(type,"two-sided")) (1-conf.level)/2 else 1-conf.level
za <- qnorm(a,lower.tail=FALSE)
a1 <- pnorm(z0+(z0-za)/(1-a0*(z0-za)))
a2 <- pnorm(z0+(z0+za)/(1-a0*(z0+za)))
q <- signif(quantile(z,p=c(a1,a2)),digits=4)
li <- c(paste("(",q[1],",",q[2],")",sep=""),paste(q[1],"(LCB)"),paste(q[2],"(UCB)"))
CI <- li[ti==type]
lims <- list(q,q[1],q[2])
llims <- lims[[c(1:3)[ti==type]]]
# FOR THE RESULTS
param <- "standard deviation"
stat.name <- "ratio.sd"
if (!test) {alt <- NULL; p.value <- NULL; p <- NULL}
results <- 
 list(Stacked=stacked,Boot.values=z,Confidence.limits=llims,Parameter=param,
      Header=paste("RESULTS OF BCa BOOTSTRAP FOR",toupper(stat.name)),
      Variable=variable,Pop.1=pop.1,Pop.2=pop.2,n.1=m,n.2=n,Statistic=stat.name,
      Observed.1=obs.1,Observed.2=obs.2,Observed=obs,Replications=R,Mean=mean(z),
      SE=sd(z),Bias=bias,Percent.bias=percent.bias,Null=null.hyp,Alternative=alt,
      P.value=p.value,p.value=p,Level=cl,Type=type,Confidence.interval=CI)
class(results) <- "boot.two"  # bootstrap, two sample.
results
}
