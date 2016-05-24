boot.one.bca <-
function(x, parameter, null.hyp = NULL,
                         alternative = c("two.sided","less","greater"),
                         conf.level = 0.95, type = NULL, R = 9999)
{
# require(boot)
proportion <- mean
test <- !is.null(null.hyp)
statistic <- parameter
boot.statistic <- boot(x,function(x,i) statistic(x[i]),R=R)
z <- c(boot.statistic$t)
bias <- signif(mean(z)-statistic(x),digits=3)
percent.bias <- signif(100*abs(bias/statistic(x)),digits=3)
cl <- paste(100*conf.level,"%",sep="")
if (identical(alternative,c("two.sided","less","greater"))) alternative <- "two.sided"
z1 <- mean(z < statistic(x))
z0 <- qnorm(z1)
n <- length(x)
d <- vector(length=n)
for (j in 1:n) d[j] <- statistic(x[-j])
mean.d <- mean(d);
a0 <- -sum((d-mean.d)^3)/(6*(sum((d-mean.d)^2))^(3/2));
# FOR THE HYPOTHESIS TEST
if (test)
{
 rtp <- (sum((z < null.hyp))+(sum((z == null.hyp))+1)/2)/(R+1)
 b0 <- qnorm(rtp)
 c0 <- ((2+a0*b0-a0*z0)*z0-b0)/(1+a0*b0-a0*z0)
 p0 <- pnorm(c0); # P-value for a left-tailed test
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
li <- c(paste("(",q[1],", ",q[2],")",sep=""),paste(q[1],"(LCB)"),paste(q[2],"(UCB)"))
CI <- li[ti==type]
lims <- list(q,q[1],q[2])
llims <- lims[[c(1:3)[ti==type]]]
# FOR THE RESULTS
var.name <- all.names(substitute(x))
if (length(var.name)>1) var.name <- var.name[[3]]
stat.name <- all.names(substitute(parameter));
if (!test) {alt <- NULL; p.value <- NULL; p <- NULL}
results <- 
 list(Boot.values=z,Confidence.limits=llims,
      Header=paste("RESULTS OF BCa BOOTSTRAP FOR",toupper(stat.name)),
      Variable=var.name,n=n,Statistic=stat.name,Observed=statistic(x),
      Replications=R,Mean=mean(z),SE=sd(z),Bias=bias,
      Percent.bias=percent.bias,Null=null.hyp,Alternative=alt,P.value=p.value,
      p.value=p,Level=cl,Type=type,Confidence.interval=CI);
class(results) <- "boot.one"  # bootstrap, one sample.
results
}
