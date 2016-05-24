`etc.diff` <-
function(formula,data,base=1,margin.up=NULL,margin.lo=-margin.up,
                     method="var.unequal",FWER=0.05) {


if (length(formula) != 3) {
  stop("formula mis-specified")
}
mf <- model.frame(formula, data)
if (ncol(mf) != 2) {
  stop("Specify one response and only one class variable in the formula")
}
if (is.numeric(mf[, 1]) == FALSE) {
  stop("Response variable must be numeric")
}
Response <- mf[, 1]
Treatment <- as.factor(mf[, 2])
tr.names <- levels(Treatment)
comp.names <- paste(tr.names[-base], tr.names[base], sep = "-")
k <- length(comp.names)                                                                         # number of comparisons

if ( is.numeric(margin.up)==FALSE | (length(margin.up)==k)+(length(margin.up)==1)==0 ) {
  stop("margin.up must be a single numeric value or a numeric vector of lenght equal to the number of comparisons")
}
if (length(margin.up)==1) {
  margin.up <- rep(margin.up,k)
}
if ( is.numeric(margin.lo)==FALSE | (length(margin.lo)==k)+(length(margin.lo)==1)==0 ) {
  stop("margin.lo must be a single numeric value or a numeric vector of lenght equal to the number of comparisons")
}
if (length(margin.lo)==1) {
  margin.lo <- rep(margin.lo,k)
}
if (any(margin.up<=0) | any(margin.lo>=0)) {
  stop("All components of margin.up (margin.lo) must be positiv (negative)")
}

method <- match.arg(method, choices = c("Bofinger", "var.equal", "var.unequal", "non.par"))

tr.mean   <- tapply(Response,Treatment,mean)
tr.sd     <- tapply(Response,Treatment,sd)
tr.n      <- tapply(Response,Treatment,length)
estimate  <- tr.mean[-base]-tr.mean[base]                                                       # estimates
test.stat <- numeric(k)
m <- floor(k/2); u <- m+1
p.value <- numeric(k)

if (method=="Bofinger")                                                                         # due to Bof./Tong, only exact
{                                                                                               # for balancedness!
  if (any(margin.up!=-margin.lo)) {
    stop("Method Bofinger works only for margin.up = -margin.lo")
  }
  if (all(as.numeric(tr.n[-base])==tr.n[-base][1])==FALSE)
  {
    cat("Warning: Method Bofinger is only correct for equal sample sizes of the test treatments", 
        "\n")
  }
  s <- sqrt( sum((tr.n-1)*tr.sd^2)/sum(tr.n-1) )                                                # pooled standard deviation
  degr.fr <- sum(tr.n-1)                                                                        # degree of freedom
  corr.mat <- diag(k)                                                                           # correl. matrix due to Bof./Tong
  if (k>1)
  {
    for(i in 1:k) { for(j in 1:k) { corr.mat[i,j]=1/sqrt( (1+tr.n[base]/tr.n[-base][i])*
                                                          (1+tr.n[base]/tr.n[-base][j]) ) }}
    for(i in 1:m) { for(j in u:k) { corr.mat[i,j]=-corr.mat[i,j] }}
    for(i in u:k) { for(j in 1:m) { corr.mat[i,j]=-corr.mat[i,j] }}
    diag(corr.mat)=rep(1,times=ncol(corr.mat))
  }
  qu <- qmvt(1-FWER, tail="lower.tail", df=degr.fr, corr=corr.mat)$quantile
  test.stat <- ( abs(tr.mean[-base]-tr.mean[base])-margin.up ) / ( s * sqrt(1/tr.n[-base] + 1/tr.n[base]) )
  for (i in 1:k) {
    p.value[i]=1-pmvt(lower=rep(test.stat[i],times=k),upper=Inf,df=degr.fr,corr=corr.mat)[1]
  }
  lower <- estimate-qu*( s * sqrt(1/tr.n[-base] + 1/tr.n[base]) ); lower[lower>0]=0
  upper <- estimate+qu*( s * sqrt(1/tr.n[-base] + 1/tr.n[base]) ); upper[upper<0]=0
  conf.int <- rbind(lower,upper); rownames(conf.int) <- c("lower","upper")
  value <- list(comp.names=comp.names,estimate=estimate,degr.fr=degr.fr,test.stat=test.stat,crit.value=-qu,corr.mat=corr.mat,
                p.value=p.value,conf.int=conf.int,base=base,margin.lo=margin.lo,margin.up=margin.up,method=method,
                FWER=FWER)
}

if(method=="var.equal")                                                                         # Bonferroni-adjustment
{
  s <- sqrt(sum((tr.n-1)*tr.sd^2)/sum(tr.n-1))                                                  # !: pooled standard deviation
  degr.fr <- sum(tr.n-1)                                                                        # degree of freedom
  qu <- qt(1-FWER/k, df=degr.fr, lower.tail=TRUE)
  test.stat.up <- ( tr.mean[-base]-tr.mean[base]-margin.lo ) / 
                  ( s * sqrt(1/tr.n[-base] + 1/tr.n[base]) )                                    # test "up"
  test.stat.do <- ( tr.mean[-base]-tr.mean[base]-margin.up ) / 
                  ( s * sqrt(1/tr.n[-base] + 1/tr.n[base]) )                                    # test "down"
  for (i in 1:k) {
    test.stat[i]=max(-test.stat.up[i],test.stat.do[i])
    p.value[i]=min(pt(q=test.stat[i], df=degr.fr, lower.tail=TRUE)*k, 1)
  }
  lower <- estimate-qu*( s * sqrt(1/tr.n[-base] + 1/tr.n[base]) ); lower[lower>0]=0
  upper <- estimate+qu*( s * sqrt(1/tr.n[-base] + 1/tr.n[base]) ); upper[upper<0]=0
  conf.int <- rbind(lower,upper); rownames(conf.int) <- c("lower","upper")
  value <- list(comp.names=comp.names,estimate=estimate,degr.fr=degr.fr,test.stat=test.stat,crit.value=-qu,
                p.value=p.value,conf.int=conf.int,base=base,margin.lo=margin.lo,margin.up=margin.up,method=method,
                FWER=FWER)
}

if(method=="var.unequal")                                                                       # Bonferroni-adjustment
{
  degr.fr <- ( (tr.sd[-base])^2/tr.n[-base]+(tr.sd[base])^2/tr.n[base] )^2 /                    # degrees of freedom (Welch)
             ( ((tr.sd[-base])^2/tr.n[-base])^2/(tr.n[-base]-1) + ((tr.sd[base])^2/tr.n[base])^2/(tr.n[base]-1) )
  test.stat.up <- ( tr.mean[-base]-tr.mean[base]-margin.lo ) /
                  ( sqrt((tr.sd[-base])^2/tr.n[-base] + (tr.sd[base])^2/tr.n[base]) )           # test "up"
  test.stat.do <- ( tr.mean[-base]-tr.mean[base]-margin.up ) /
                  ( sqrt((tr.sd[-base])^2/tr.n[-base] + (tr.sd[base])^2/tr.n[base]) )           # test "down"
  qu <- numeric(k)
  for (i in 1:k) {
    qu[i]=qt(1-FWER/k, df=degr.fr[i], lower.tail=TRUE)
    test.stat[i]=max(-test.stat.up[i],test.stat.do[i])
    p.value[i]=min(pt(q=test.stat[i], df=degr.fr[i], lower.tail=TRUE)*k, 1)
  }
  lower <- estimate-qu*( sqrt((tr.sd[-base])^2/tr.n[-base] + (tr.sd[base])^2/tr.n[base]) ); lower[lower>0]=0
  upper <- estimate+qu*( sqrt((tr.sd[-base])^2/tr.n[-base] + (tr.sd[base])^2/tr.n[base]) ); upper[upper<0]=0
  conf.int <- rbind(lower,upper); rownames(conf.int) <- c("lower","upper")
  value <- list(comp.names=comp.names,estimate=estimate,degr.fr=degr.fr,test.stat=test.stat,crit.value=-qu,
                p.value=p.value,conf.int=conf.int,base=base,margin.lo=margin.lo,margin.up=margin.up,method=method,
                FWER=FWER)
}

if(method=="non.par")                                                                           # Bonferroni-adjustment
{
  test.stat.up <- p.value.up <- lower <- numeric(k)
  test.stat.do <- p.value.do <- upper <- numeric(k)
  for (i in 1:k) {
    test.up <- wilcox.test(x=subset(mf,mf[,2]==tr.names[-base][i])[,1],y=subset(mf,mf[,2]==tr.names[base])[,1],
               alternative="greater",mu=margin.lo[i],paired=FALSE,exact=FALSE,correct=TRUE,conf.int=TRUE,conf.level=1-FWER/k)
    test.do <- wilcox.test(x=subset(mf,mf[,2]==tr.names[-base][i])[,1],y=subset(mf,mf[,2]==tr.names[base])[,1],
               alternative="less",mu=margin.up[i],paired=FALSE,exact=FALSE,correct=TRUE,conf.int=TRUE,conf.level=1-FWER/k)
    test.stat.up[i]=test.up$statistic; test.stat.do[i]=test.do$statistic
    p.value.up[i]=test.up$p.value; p.value.do[i]=test.do$p.value
    p.value[i]=min(max(p.value.up[i],p.value.do[i])*k,1)
    lower[i]=test.up$conf.int[1]; upper[i]=test.do$conf.int[2]
  }
  test.stat <- cbind(test.stat.up,test.stat.do)
  lower[lower>0]=0; upper[upper<0]=0
  conf.int <- rbind(lower,upper); rownames(conf.int) <- c("lower","upper")
  value <- list(comp.names=comp.names,estimate=estimate,test.stat=test.stat,
                p.value=p.value,conf.int=conf.int,base=base,margin.lo=margin.lo,margin.up=margin.up,method=method,
                FWER=FWER)
}

if(method=="var.unequal") {
  names(value$degr.fr) <- comp.names
  names(value$crit.value) <- comp.names
}
names(value$estimate) <- comp.names
names(value$test.stat) <- comp.names
names(value$p.value) <- comp.names
colnames(value$conf.int) <- comp.names
class(value) <- "etc.diff"

return(value)


}

