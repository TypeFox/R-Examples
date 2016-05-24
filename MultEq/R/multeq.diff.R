multeq.diff <-
function(data,grp,resp=NULL,base=1,margin.lo=NULL,margin.up=NULL,
                        method="single.step",var.equal=FALSE,FWER=0.05) {


if (length(grp) > 1) {
  stop("Specify only one grouping variable")
}
tr.names <- levels(data[,grp])
if (length(tr.names) > 2) {
  stop("Grouping factor must have exactly 2 levels")
}
comp.name <- paste(tr.names[-base], tr.names[base], sep = "-")
Resp.X <- subset(data, data[,grp]==tr.names[-base])                      # treatment/minuend
Resp.Y <- subset(data, data[,grp]==tr.names[base])                       # base/control/subtrahend
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
test.stat <- numeric(n.ep)

X.mean <- colMeans(Resp.X); Y.mean <- colMeans(Resp.Y)                   # mean vectors for X and Y
estimate <- X.mean-Y.mean
cov.matX <- cov(Resp.X)                                                  # just the variances needed
cov.matY <- cov(Resp.Y)                                                  # just the variances needed
if (var.equal==TRUE) {
  cov.mat=((X.n-1)*cov.matX+(Y.n-1)*cov.matY)/(X.n+Y.n-2)                # common estimated covariance matrix of the data
  degr.fr <- X.n+Y.n-2                                                   # df (single number)
} else {
  degr.fr <- ( diag(cov.matX)/X.n + diag(cov.matY)/Y.n )^2 /             # df (vector)
             ( (diag(cov.matX)/X.n)^2 / (X.n-1) + (diag(cov.matY)/Y.n)^2 / (Y.n-1) )
}

if (is.numeric(margin.lo)) {
  if (var.equal==TRUE) T.up <- (X.mean-Y.mean-margin.lo)/(  sqrt(diag(cov.mat))*sqrt( 1/X.n + 1/Y.n )  )
  else T.up <- (X.mean-Y.mean-margin.lo)/(  sqrt( diag(cov.matX)/X.n + diag(cov.matY)/Y.n )  )
}
if (is.numeric(margin.up)) {
  if (var.equal==TRUE) T.do <- (X.mean-Y.mean-margin.up)/(  sqrt(diag(cov.mat))*sqrt( 1/X.n + 1/Y.n )  )
  else T.do <- (X.mean-Y.mean-margin.up)/(  sqrt( diag(cov.matX)/X.n + diag(cov.matY)/Y.n )  )
}

if ((is.numeric(margin.lo)) & is.numeric(margin.up)) {
  for (i in 1:n.ep) test.stat[i] <- min(T.up[i],-T.do[i])
}
if ((is.numeric(margin.lo)) & is.numeric(margin.up)==FALSE) {
  test.stat <- T.up
}
if ((is.numeric(margin.lo)==FALSE) & is.numeric(margin.up)) {
  test.stat <- -T.do
}

if (method == "step.up") {
  p.value <- numeric(n.ep)
  lower <- numeric(n.ep); upper <- numeric(n.ep)
  if (var.equal==TRUE) {
    for (i in 1:n.ep) {
      pos <- which(p.value<FWER)
      if (length(pos)==n.ep+1-i){
        p.value[pos]=pt(q=test.stat[pos],df=degr.fr,lower.tail=FALSE)*i
        if (is.numeric(margin.lo)) {
          lower[pos]=X.mean[pos]-Y.mean[pos]-qt(p=1-FWER/i,df=degr.fr)*sqrt(diag(cov.mat)[pos])*sqrt( 1/X.n + 1/Y.n )
        } else { lower[pos]=-Inf }
        if (is.numeric(margin.up)) {
          upper[pos]=X.mean[pos]-Y.mean[pos]+qt(p=1-FWER/i,df=degr.fr)*sqrt(diag(cov.mat)[pos])*sqrt( 1/X.n + 1/Y.n )
        } else { upper[pos]=Inf }
      }
    }
  } else {
    for (i in 1:n.ep) {
      pos <- which(p.value<FWER)
      if (length(pos)==n.ep+1-i){
        p.value[pos]=pt(q=test.stat[pos],df=degr.fr[pos],lower.tail=FALSE)*i
        if (is.numeric(margin.lo)) {
          lower[pos]=X.mean[pos]-Y.mean[pos]-qt(p=1-FWER/i,df=degr.fr[pos])*
                     sqrt( diag(cov.matX)[pos]/X.n + diag(cov.matY)[pos]/Y.n )
        } else { lower[pos]=-Inf }
        if (is.numeric(margin.up)) {
          upper[pos]=X.mean[pos]-Y.mean[pos]+qt(p=1-FWER/i,df=degr.fr[pos])*
                     sqrt( diag(cov.matX)[pos]/X.n + diag(cov.matY)[pos]/Y.n )
        } else { upper[pos]=Inf }
      }
    }
  }
}
if (method == "single.step")
{
  p.value=pt(q=test.stat,df=degr.fr,lower.tail=FALSE)*n.ep
  if (var.equal==TRUE) {
    if (is.numeric(margin.lo)) {
      lower=X.mean-Y.mean-qt(p=1-FWER/n.ep,df=degr.fr)*sqrt(diag(cov.mat))*sqrt( 1/X.n + 1/Y.n )
    } else { lower=rep(-Inf,n.ep) }
    if (is.numeric(margin.up)) {
      upper=X.mean-Y.mean+qt(p=1-FWER/n.ep,df=degr.fr)*sqrt(diag(cov.mat))*sqrt( 1/X.n + 1/Y.n )
    } else { upper=rep(Inf,n.ep) }
  } else {
    if (is.numeric(margin.lo)) {
      lower=X.mean-Y.mean-qt(p=1-FWER/n.ep,df=degr.fr)*sqrt( diag(cov.matX)/X.n + diag(cov.matY)/Y.n )
    } else { lower=rep(-Inf,n.ep) }
    if (is.numeric(margin.up)) {
      upper=X.mean-Y.mean+qt(p=1-FWER/n.ep,df=degr.fr)*sqrt( diag(cov.matX)/X.n + diag(cov.matY)/Y.n )
    } else { upper=rep(Inf,n.ep) }
  }
}
p.value[p.value>1]=1

value <- list(comp.name=comp.name,estimate=estimate,degr.fr=degr.fr,test.stat=test.stat,
              p.value=p.value,lower=lower,upper=upper,margin.lo=margin.lo,margin.up=margin.up,
              base=base,method=method,var.equal=var.equal,FWER=FWER)

names(value$test.stat) <- names(estimate)
names(value$p.value) <- names(estimate)
names(value$lower) <- names(value$upper) <- names(estimate)
if (is.numeric(margin.lo)) { names(value$margin.lo) <- names(estimate) }
if (is.numeric(margin.up)) { names(value$margin.up) <- names(estimate) }
class(value) <- "multeq.diff"

return(value)


}

