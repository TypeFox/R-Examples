boot.cond.mean.per <-
function(x, y, xp, conf.level = 0.95, R = 9999)
{
# require(boot)
test <- FALSE
null.hyp <- NULL; type <- NULL
data <- data.frame(x,y)
dp <- data.frame(x=xp)
obs <- predict(lm(y~x),dp)
boot.fun <- function(d,i) {u <- d$x[i]; v <- d$y[i]; predict(lm(v~u),data.frame(u=xp))}
boot.predict <- suppressWarnings(boot(data,boot.fun,R=R))
z <- c(boot.predict$t)
z <- z[(!is.infinite(z) & !is.na(z))] # To eliminate resamples with undefined slope due to a 0 variance.
R <- length(z)  # To adjust for eliminated resamples.
bias <- signif(mean(z)-obs,digits=3)
percent.bias <- signif(100*abs(bias/obs),digits=3)
cl <- paste(100*conf.level,"%",sep="")
# FOR THE CONFIDENCE INTERVAL
a1 <- (1-conf.level)/2 
a2 <- 1-a1
q <- signif(quantile(z,p=c(a1,a2)),digits=4)
CI <- paste("(",q[1],", ",q[2],")",sep="")
# FOR THE RESULTS
var1.name <- all.names(substitute(x))
if (length(var1.name)>1) var1.name <- var1.name[[3]]
var2.name <- all.names(substitute(y))
if (length(var2.name)>1) var2.name <- var2.name[[3]]
stat.name <- "conditional mean"
alt <- NULL; p.value <- NULL; p <- NULL
results <- 
 list(Boot.values=z,Confidence.limits=q,
      Header=paste("RESULTS OF PERCENTILE BOOTSTRAP FOR",toupper(stat.name),"CORRESPONDING TO",var1.name,"=",xp),
      Variable.1=var1.name,Variable.2=var2.name,n=length(x),
      Statistic="Fit",Observed=obs,Replications=R,Mean=mean(z),
      SE=sd(z),Bias=bias,Percent.bias=percent.bias,Null=null.hyp,
      Alternative=alt,P.value=p.value,p.value=p,Level=cl,Type=type,
      Confidence.interval=CI,cor.ana=FALSE)
class(results) <- "boot.regcor"  # bootstrap, regression and correlation
results
}
