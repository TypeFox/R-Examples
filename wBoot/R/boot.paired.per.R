boot.paired.per <-
function(x, y, variable=NULL, null.hyp = NULL,
                            alternative = c("two.sided","less","greater"), conf.level = 0.95,
                            type = NULL, R = 9999)
{
# require(boot)
test <- !is.null(null.hyp)
pop.1 <- all.names(substitute(x));
if (length(pop.1)>1) pop.1 <- paste(pop.1[[2]],pop.1[[1]],pop.1[[3]],sep="")
pop.2 <- all.names(substitute(y));
if (length(pop.2)>1) pop.2 <- paste(pop.2[[2]],pop.2[[1]],pop.2[[3]],sep="")
data <- x-y
obs <- mean(data)
boot.mean <- boot(data,function(w,i) mean(w[i]),R=R)
z <- c(boot.mean$t)
bias <- signif(mean(z)-obs,digits=3)
percent.bias <- signif(100*abs(bias/obs),digits=3)
cl <- paste(100*conf.level,"%",sep="")
if (identical(alternative,c("two.sided","less","greater"))) alternative <- "two.sided"
# FOR THE HYPOTHESIS TEST
if (test)
{
 p0 <- (sum((z > null.hyp))+(sum((z == null.hyp))+1)/2)/(R+1) # P-value for a left-tailed test
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
a1 <- if (identical(type,"two-sided")) (1-conf.level)/2 else 1-conf.level
a2 <- 1-a1
q <- signif(quantile(z,p=c(a1,a2)),digits=4)
li <- c(paste("(",q[1],", ",q[2],")",sep=""),paste(q[1],"(LCB)"),paste(q[2],"(UCB)"))
CI <- li[ti==type]
lims <- list(q,q[1],q[2])
llims <- lims[[c(1:3)[ti==type]]]
# FOR THE RESULTS
if (!test) {alt <- NULL; p.value <- NULL; p <- NULL}
stat.name <- "diff.mean"
results <- 
 list(Boot.values=z,Confidence.limits=llims,
      Header=paste("RESULTS OF PERCENTILE BOOTSTRAP FOR",toupper(stat.name)),
      Variable=variable,Pop.1=pop.1,Pop.2=pop.2,n=length(data),
      Statistic=stat.name,Observed=obs,Replications=R,Mean=mean(z),
      SE=sd(z),Bias=bias,Percent.bias=percent.bias,Null=null.hyp,
      Alternative=alt,P.value=p.value,p.value=p,Level=cl,Type=type,
      Confidence.interval=CI)
class(results) <- "boot.paired"  # bootstrap, paired sample.
results
}
