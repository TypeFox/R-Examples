multeq.rat <-
function(data,grp,resp=NULL,base=1,margin.lo=NULL,margin.up=NULL,
                       method="single.step",var.equal=FALSE,FWER=0.05) {


if (length(grp) > 1) {
  stop("Specify only one grouping variable")
}
tr.names <- levels(data[,grp])
if (length(tr.names) > 2) {
  stop("Grouping factor must have exactly 2 levels")
}
comp.name <- paste(tr.names[-base], tr.names[base], sep = "/")
Resp.X <- subset(data, data[,grp]==tr.names[-base])                      # treatment/numerator
Resp.Y <- subset(data, data[,grp]==tr.names[base])                       # base/control/denominator
if (is.null(resp)) {
  n.ep <- length(names(data))-1                                          # number of endpoints
  Resp.X <- Resp.X[,-which(names(data)%in%grp)]
  Resp.Y <- Resp.Y[,-which(names(data)%in%grp)]
} else {
  n.ep <- length(resp)                                                   # number of endpoints
  Resp.X <- Resp.X[resp]
  Resp.Y <- Resp.Y[resp]
}

if (is.numeric(margin.lo) & length(margin.lo) != n.ep) {
  stop("Length of margin.lo is not equal to the number of response variables")
}
if (is.numeric(margin.up) & length(margin.up) != n.ep) {
  stop("Length of margin.up is not equal to the number of response variables")
}
method <- match.arg(method, choices = c("single.step", "step.up"))

if (!is.logical(var.equal)) { 
  stop("var.equal must be TRUE or FALSE")
}

X.n <- nrow(Resp.X); Y.n <- nrow(Resp.Y)                                 # sample sizes for X and Y
X.mean <- colMeans(Resp.X); Y.mean <- colMeans(Resp.Y)                   # mean vectors for X and Y
estimate <- X.mean/Y.mean
cov.matX <- cov(Resp.X)                                                  # just the variances needed
cov.matY <- cov(Resp.Y)                                                  # just the variances needed
if (var.equal==TRUE) {
  cov.mat=((X.n-1)*cov.matX+(Y.n-1)*cov.matY)/(X.n+Y.n-2)                # common estimated covariance matrix of the data
  degr.fr <- X.n+Y.n-2                                                   # df (single number)
} else {
  degr.fr.up <- ( diag(cov.matX)/X.n + margin.lo^2*diag(cov.matY)/Y.n )^2 /   # Welch degrees of freedom for test "up"
                ( (diag(cov.matX)/X.n)^2 / (X.n-1) + (margin.lo^2*diag(cov.matY)/Y.n)^2 / (Y.n-1) )
  degr.fr.do <- ( diag(cov.matX)/X.n + margin.up^2*diag(cov.matY)/Y.n )^2 /   # Welch degrees of freedom for test "do"
                ( (diag(cov.matX)/X.n)^2 / (X.n-1) + (margin.up^2*diag(cov.matY)/Y.n)^2 / (Y.n-1) )
  degr.fr.ci <- ( diag(cov.matX)/X.n + estimate^2*diag(cov.matY)/Y.n )^2 /    # Welch degrees of freedom for CI (plug in)
                ( (diag(cov.matX)/X.n)^2 / (X.n-1) + (estimate^2*diag(cov.matY)/Y.n)^2 / (Y.n-1) )
}

if (is.numeric(margin.lo)) {
  if (var.equal==TRUE) test.stat.up <- (X.mean-margin.lo*Y.mean)/(  sqrt(diag(cov.mat))*sqrt( 1/X.n + margin.lo^2/Y.n )  )
  else  test.stat.up <- (X.mean-margin.lo*Y.mean)/(  sqrt( diag(cov.matX)/X.n + diag(cov.matY)*margin.lo^2/Y.n )  )
}
if (is.numeric(margin.up)) {
  if (var.equal==TRUE) test.stat.do <- (X.mean-margin.up*Y.mean)/(  sqrt(diag(cov.mat))*sqrt( 1/X.n + margin.up^2/Y.n )  )
  else test.stat.do <- (X.mean-margin.up*Y.mean)/(  sqrt( diag(cov.matX)/X.n + diag(cov.matY)*margin.up^2/Y.n )  )
}

if (var.equal==TRUE)
{
  test.stat <- numeric(n.ep)
  if ((is.numeric(margin.lo)) & is.numeric(margin.up)) {
    for (i in 1:n.ep) test.stat[i] <- min(test.stat.up[i],-test.stat.do[i])
  }
  if ((is.numeric(margin.lo)) & is.numeric(margin.up)==FALSE) {
    test.stat <- test.stat.up
  }
  if ((is.numeric(margin.lo)==FALSE) & is.numeric(margin.up)) {
    test.stat <- -test.stat.do
  }
}

p.value <- numeric(n.ep)

if (method == "step.up") {
  lower <- numeric(n.ep); upper <- numeric(n.ep)
  if (var.equal==TRUE) {
    g <- numeric(n.ep)                                                   # see Dilba et al. 2004 p. 446
    for (i in 1:n.ep){
      pos <- which(p.value<FWER)
      if (length(pos)==n.ep+1-i){
        p.value[pos]=pt(q=test.stat[pos],df=degr.fr,lower.tail=FALSE)*i
        g[pos]=(qt(p=1-FWER/i,df=degr.fr))^2*diag(cov.mat)[pos]/( Y.n*(Y.mean[pos])^2 )
        if ( all(g<1) ) {
          if (is.numeric(margin.lo)) {
            lower[pos]=( estimate[pos]-sqrt(g[pos])*sqrt( (estimate[pos])^2+(1-g[pos])*Y.n/X.n ) )/(1-g[pos])
          } else { lower[pos]=-Inf }
          if (is.numeric(margin.up)) {
            upper[pos]=( estimate[pos]+sqrt(g[pos])*sqrt( (estimate[pos])^2+(1-g[pos])*Y.n/X.n ) )/(1-g[pos])
          } else { upper[pos]=Inf }
        } else { lower <- upper <- "NSD" }
      }
    }
  } else {
    p.value.up <- numeric(n.ep); p.value.do <- numeric(n.ep)
    for (i in 1:n.ep) {
      pos <- which(p.value<FWER)
      if (length(pos)==n.ep+1-i){
        p.value.up[pos]=pt(q=test.stat.up[pos],df=degr.fr.up[pos],lower.tail=FALSE)*i
        p.value.do[pos]=pt(q=test.stat.do[pos],df=degr.fr.do[pos],lower.tail=TRUE)*i
        for (i in seq(along.with=pos)) p.value[i] <- max(p.value.up[i],p.value.do[i])
        Ai <- (Y.mean[pos])^2-qt(p=1-FWER/i,df=degr.fr.ci[pos])^2*diag(cov.matY)[pos]/Y.n
        Bi <- -2*X.mean[pos]*Y.mean[pos]
        Ci <- (X.mean[pos])^2-qt(p=1-FWER/i,df=degr.fr.ci[pos])^2*diag(cov.matX)[pos]/X.n
        Discrimi <- Bi^2-4*Ai*Ci
        if ( all(Ai > 0) & all(Discrimi >= 0) ) {
          if (is.numeric(margin.lo)) {
            lower[pos]=( -Bi-sqrt(Discrimi) ) / ( 2*Ai )
          } else { lower[pos]=-Inf }
          if (is.numeric(margin.up)) {
            upper[pos]=( -Bi+sqrt(Discrimi) ) / ( 2*Ai )
          } else { upper[pos]=Inf }
        } else lower <- upper <- "NSD"
      }
    }
  }
}
if (method == "single.step")
{
  if (var.equal==TRUE) {
    p.value <- pt(q=test.stat,df=degr.fr,lower.tail=FALSE)*n.ep
    g=(qt(p=1-FWER/n.ep,df=degr.fr))^2*diag(cov.mat)/( Y.n*(Y.mean)^2 )
    if ( all(g<1) ) {
      if (is.numeric(margin.lo)) {
        lower=( estimate-sqrt(g)*sqrt( (estimate)^2+(1-g)*Y.n/X.n ) )/(1-g)
      } else { lower=rep(-Inf,n.ep) }
      if (is.numeric(margin.up)) {
        upper=( estimate+sqrt(g)*sqrt( (estimate)^2+(1-g)*Y.n/X.n ) )/(1-g)
      } else { upper=rep(Inf,n.ep) }
    } else { lower <- upper <- "NSD" }
  } else {
    p.value.up <- pt(q=test.stat.up,df=degr.fr.up,lower.tail=FALSE)*n.ep
    p.value.do <- pt(q=test.stat.do,df=degr.fr.do,lower.tail=TRUE)*n.ep
    for (i in 1:n.ep) p.value[i] <- max(p.value.up[i],p.value.do[i])
    Ai <- (Y.mean)^2-qt(p=1-FWER/n.ep,df=degr.fr.ci)^2*diag(cov.matY)/Y.n
    Bi <- -2*X.mean*Y.mean
    Ci <- (X.mean)^2-qt(p=1-FWER/n.ep,df=degr.fr.ci)^2*diag(cov.matX)/X.n
    Discrimi <- Bi^2-4*Ai*Ci
    if ( all(Ai > 0) & all(Discrimi >= 0) ) {
      if (is.numeric(margin.lo)) {
        lower=( -Bi-sqrt(Discrimi) ) / ( 2*Ai )
      } else { lower=rep(-Inf,n.ep) }
      if (is.numeric(margin.up)) {
        upper=( -Bi+sqrt(Discrimi) ) / ( 2*Ai )
      } else { upper=rep(Inf,n.ep) }
    } else lower <- upper <- "NSD"
  }
}
p.value[p.value>1]=1

if (var.equal==TRUE)
{
  value <- list(comp.name=comp.name,estimate=estimate,degr.fr=degr.fr,test.stat=test.stat,
                p.value=p.value,lower=lower,upper=upper,margin.lo=margin.lo,margin.up=margin.up,
                base=base,method=method,var.equal=var.equal,FWER=FWER)
  names(value$test.stat) <- names(estimate)
} else {
  value <- list(comp.name=comp.name,estimate=estimate,degr.fr.up=degr.fr.up,degr.fr.do=degr.fr.do,
                degr.fr.ci=degr.fr.ci,test.stat.up=test.stat.up,test.stat.do=test.stat.do,p.value=p.value,
                lower=lower,upper=upper,margin.lo=margin.lo,margin.up=margin.up,
                base=base,method=method,var.equal=var.equal,FWER=FWER)
  names(value$test.stat.up) <- names(value$test.stat.do) <- names(estimate)
}

names(value$p.value) <- names(estimate)
if (is.numeric(lower)) { names(value$lower) <- names(value$upper) <- names(estimate) } # both together
if (is.numeric(margin.lo)) { names(value$margin.lo) <- names(estimate) }
if (is.numeric(margin.up)) { names(value$margin.up) <- names(estimate) }
class(value) <- "multeq.rat"

return(value)


}

