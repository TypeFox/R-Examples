boot.cor.per <-
function(x, y, null.hyp = NULL, alternative = c("two.sided","less","greater"),
                         conf.level = 0.95, type = NULL, R = 9999)
{
# require(boot)
obs <- cor(x,y)
test <- !is.null(null.hyp)
data <- data.frame(x,y)
boot.cor <- suppressWarnings(boot(data,function(d,i) cor(d$x[i],d$y[i]),R=R))
z <- c(boot.cor$t)
z <- z[(!is.infinite(z) & !is.na(z))] # To eliminate resamples with undefined correlation due to a 0 variance.
R <- length(z)  # To adjust for eliminated resamples.
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
var1.name <- all.names(substitute(x))
if (length(var1.name)>1) var1.name <- var1.name[[3]]
var2.name <- all.names(substitute(y))
if (length(var2.name)>1) var2.name <- var2.name[[3]]
stat.name <- "correlation"
if (!test) {alt <- NULL; p.value <- NULL; p <- NULL}
results <- 
 list(Boot.values=z,Confidence.limits=llims,
      Header=paste("RESULTS OF PERCENTILE BOOTSTRAP FOR",toupper(stat.name)),
      Variable.1=var1.name,Variable.2=var2.name,n=length(x),
      Statistic=stat.name,Observed=obs,Replications=R,Mean=mean(z),
      SE=sd(z),Bias=bias,Percent.bias=percent.bias,Null=null.hyp,
      Alternative=alt,P.value=p.value,p.value=p,Level=cl,Type=type,
      Confidence.interval=CI,cor.ana=TRUE)
class(results) <- "boot.regcor"  # bootstrap, regression and correlation
results
}
