`etc.rat` <-
function(formula,data,base=1,margin.up=NULL,margin.lo=1/margin.up,
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
comp.names <- paste(tr.names[-base], tr.names[base], sep = "/")
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
if (any(margin.up<=1) | any(margin.lo>=1)) {
  stop("All components of margin.up (margin.lo) must be larger (smaller) than 1")
}
method <- match.arg(method, choices = c("var.equal", "var.unequal"))

tr.mean  <- tapply(Response,Treatment,mean)
tr.sd    <- tapply(Response,Treatment,sd)
tr.n     <- tapply(Response,Treatment,length)
estimate <- tr.mean[-base]/tr.mean[base]                                                        # estimates
m <- floor(k/2); u <- m+1
test.stat <- numeric(k)
p.value <- numeric(k)

if(method=="var.equal")                                                                         # Bonferroni-adjustment
{
  s <- sqrt(sum((tr.n-1)*tr.sd^2)/sum(tr.n-1))                                                  # !: pooled standard deviation
  degr.fr <- sum(tr.n-1)                                                                        # degree of freedom
  qu <- qt(1-FWER/k, df=degr.fr, lower.tail=TRUE)
  test.stat.up <- ( tr.mean[-base]-margin.lo*tr.mean[base] ) /
                  ( s * sqrt(1/tr.n[-base] + margin.lo^2/tr.n[base]) )                          # test "up"
  test.stat.do <- ( tr.mean[-base]-margin.up*tr.mean[base] ) /
                  ( s * sqrt(1/tr.n[-base] + margin.up^2/tr.n[base]) )                          # test "down"
  for (i in 1:k) {
    test.stat[i]=max(-test.stat.up[i],test.stat.do[i])
    p.value[i]=min(pt(q=test.stat[i], df=degr.fr, lower.tail=TRUE)*k, 1)
  }
  Ai <- (tr.mean[base])^2-qu^2*s^2/tr.n[base]
  Bi <- -2*tr.mean[-base]*tr.mean[base]
  Ci <- (tr.mean[-base])^2-qu^2*s^2/tr.n[-base]
  Discrimi <- Bi^2-4*Ai*Ci
  if ( all(Ai > 0) & all(Discrimi >= 0) ) {
    lower <- ( -Bi-sqrt(Discrimi) ) / ( 2*Ai ); lower[lower>1]=1
    upper <- ( -Bi+sqrt(Discrimi) ) / ( 2*Ai ); upper[upper<1]=1
    conf.int <- rbind(lower,upper); rownames(conf.int) <- c("lower","upper")
  } else conf.int <- "NSD"
  value <- list(comp.names=comp.names,estimate=estimate,degr.fr=degr.fr,test.stat=test.stat,crit.value=-qu,
                p.value=p.value,conf.int=conf.int,base=base,margin.lo=margin.lo,margin.up=margin.up,
                method=method,FWER=FWER)
  names(value$test.stat) <- comp.names
}

if(method=="var.unequal")                                                                       # Bonferroni-adjustment
{
  degr.fr.up <- ( (tr.sd[-base])^2/tr.n[-base]+(margin.lo*tr.sd[base])^2/tr.n[base] )^2 /       # Welch degrees of freedom
                ( ((tr.sd[-base])^2/tr.n[-base])^2/(tr.n[-base]-1) +                            # for test "up"
                  ((margin.lo*tr.sd[base])^2/tr.n[base])^2/(tr.n[base]-1) )
  degr.fr.do <- ( (tr.sd[-base])^2/tr.n[-base]+(margin.up*tr.sd[base])^2/tr.n[base] )^2 /       # Welch degrees of freedom
                ( ((tr.sd[-base])^2/tr.n[-base])^2/(tr.n[-base]-1) +                            # for test "do"
                  ((margin.up*tr.sd[base])^2/tr.n[base])^2/(tr.n[base]-1) )
  degr.fr.ci <- ( (tr.sd[-base])^2/tr.n[-base]+(estimate*tr.sd[base])^2/tr.n[base] )^2 /        # Welch degrees of freedom
                ( ((tr.sd[-base])^2/tr.n[-base])^2/(tr.n[-base]-1) +                            # for CI (plug in)
                  ((estimate*tr.sd[base])^2/tr.n[base])^2/(tr.n[base]-1) )
  qu.up <- numeric(k); qu.do <- numeric(k); qu.ci <- numeric(k)
  test.stat.up <- ( tr.mean[-base]-margin.lo*tr.mean[base] ) /
                  ( sqrt(tr.sd[-base]^2/tr.n[-base] + margin.lo^2*tr.sd[base]^2/tr.n[base]) )   # test "up"
  test.stat.do <- ( tr.mean[-base]-margin.up*tr.mean[base] ) /
                  ( sqrt(tr.sd[-base]^2/tr.n[-base] + margin.up^2*tr.sd[base]^2/tr.n[base]) )   # test "down"
  p.value.up <- numeric(k); p.value.do <- numeric(k)
  for (i in 1:k) {
    qu.up[i]=qt(1-FWER/k, df=degr.fr.up[i], lower.tail=TRUE)
    qu.do[i]=qt(1-FWER/k, df=degr.fr.do[i], lower.tail=FALSE)
    qu.ci[i]=qt(1-FWER/k, df=degr.fr.ci[i], lower.tail=TRUE)
    p.value.up[i]=min(pt(q=test.stat.up[i], df=degr.fr.up[i], lower.tail=FALSE)*k, 1)
    p.value.do[i]=min(pt(q=test.stat.do[i], df=degr.fr.do[i], lower.tail=TRUE)*k, 1)
    p.value[i] <- max(p.value.up[i],p.value.do[i])
  }
  Ai <- (tr.mean[base])^2-qu.ci^2*tr.sd[base]^2/tr.n[base]
  Bi <- -2*tr.mean[-base]*tr.mean[base]
  Ci <- (tr.mean[-base])^2-qu.ci^2*tr.sd[-base]^2/tr.n[-base]
  Discrimi <- Bi^2-4*Ai*Ci
  if ( all(Ai > 0) & all(Discrimi >= 0) ) {
    lower <- ( -Bi-sqrt(Discrimi) ) / ( 2*Ai ); lower[lower>1]=1
    upper <- ( -Bi+sqrt(Discrimi) ) / ( 2*Ai ); upper[upper<1]=1
    conf.int <- rbind(lower,upper); rownames(conf.int) <- c("lower","upper")
  } else conf.int <- "NSD"
  value <- list(comp.names=comp.names,estimate=estimate,degr.fr.up=degr.fr.up,degr.fr.do=degr.fr.do,
                degr.fr.ci=degr.fr.ci,test.stat.up=test.stat.up,test.stat.do=test.stat.do,crit.value.up=qu.up,
                crit.value.do=qu.do,crit.value.ci=qu.ci,p.value=p.value,conf.int=conf.int,base=base,
                margin.lo=margin.lo,margin.up=margin.up,method=method,FWER=FWER)
  names(value$degr.fr.up) <- comp.names;    names(value$degr.fr.do) <- comp.names;    names(value$degr.fr.ci) <- comp.names
  names(value$crit.value.up) <- comp.names; names(value$crit.value.do) <- comp.names; names(value$crit.value.ci) <- comp.names
  names(value$test.stat.up) <- comp.names;  names(value$test.stat.do) <- comp.names
}

names(value$estimate) <- comp.names
names(value$p.value) <- comp.names
if (is.numeric(value$conf.int)) colnames(value$conf.int) <- comp.names
names(value$margin.up) <- comp.names; names(value$margin.lo) <- comp.names
class(value) <- "etc.rat"

return(value)


}

