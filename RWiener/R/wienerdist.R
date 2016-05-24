check_wiener_pars <- function(alpha,tau,beta,delta)
{
  if(!is.numeric(alpha) || !is.numeric(tau) || 
     !is.numeric(beta) || !is.numeric(delta)) {
    return(FALSE)
  }
  if(alpha > 0 & 
     tau > 0 &
     beta >= 0 & beta <= 1) return(TRUE)
  else return(FALSE)
}

dwiener <- function(q, alpha,tau,beta,delta, resp="upper", give_log=FALSE) 
{
  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.numeric(q) || !(is.character(resp) || is.factor(resp))) {
    stop("bad parameter values!")
  }

  if (!(length(resp) == length(q))) {
    stop("argument q and resp need to be of the same length!")
  }

  if(class(resp) == "factor") {
    resp <- as.character(resp)
  }

  d <- vector("double", length=length(q))
  for (i in 1:length(q)) {
    if (q[i]<0) stop("q must be > 0!")

    if (resp[i] == "upper") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta, give_log)
    else if (resp[i] == "lower") 
      d[i] <- .Call(dwiener_c, -q[i], alpha,tau,beta,delta, give_log)
    else if (resp[i] == "both") 
      d[i] <- .Call(dwiener_c, q[i], alpha,tau,beta,delta, give_log) +
           .Call(dwiener_c, -q[i], alpha,tau,beta,delta, give_log)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(d[i])) d[i] <- 0
  }

  return(d)
}

pwiener <- function(q, alpha,tau,beta,delta, resp="upper")
{
  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.numeric(q) || !(is.character(resp) || is.factor(resp))) {
    stop("bad parameter values!")
  }

  if (!(length(resp) == length(q))) {
    stop("argument q and resp need to be of the same length!")
  }

  if(class(resp) == "factor") {
    resp <- as.character(resp)
  }

  p <- vector("double", length=length(q))
  for (i in 1:length(q)) {
    if (q[i]<0) stop("q must be > 0!")

    if (resp[i] == "upper") 
      p[i] <- .Call(pwiener_c, q[i], alpha,tau,beta,delta)
    else if (resp[i] == "lower")
      p[i] <- .Call(pwiener_c, -q[i], alpha,tau,beta,delta)
    else if (resp[i] == "both")
      p[i] <- .Call(pwiener_full_c, q[i], alpha,tau,beta,delta)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(p[i])) p[i] <- 0
  }

  return(p)
}

qwiener <- function(p, alpha,tau,beta,delta, resp="upper")
{
  if (!check_wiener_pars(alpha,tau,beta,delta) ||
      !is.numeric(p) || !(is.character(resp) || is.factor(resp))) {
    stop("bad parameter values!")
  }

  if (!(length(resp) == length(p))) {
    stop("argument p and resp need to be of the same length!")
  }

  if(class(resp) == "factor") {
    resp <- as.character(resp)
  }


  q <- vector("double", length=length(q))
  for (i in 1:length(p)) {
    if (p[i]<0) stop("p must be > 0!")

    if (resp[i] == "upper")
      q[i] <- .Call(qwiener_c, p[i], alpha,tau,beta,delta)
    else if (resp[i] == "lower")
      q[i] <- .Call(qwiener_c, -p[i], alpha,tau,beta,delta)
    else if (resp[i] == "both")
      q[i] <- .Call(qwiener_full_c, p[i], alpha,tau,beta,delta)
    else stop("resp must be either 'lower', 'upper' or 'both'")
    if(is.nan(q[i])) p[i] <- 0
  }

  return(q)
}

rwiener <- function(n, alpha,tau,beta,delta)
{
  if (!check_wiener_pars(alpha,tau,beta,delta)) {
    stop("bad parameter values!")
  }

  rdat <- data.frame(q=vector("double"),resp=factor(levels=c("upper", "lower")))

  for (i in 1:n) {
    r <- .Call(rwiener_c, alpha,tau,beta,delta)
    if (r >= 0) rdat[i,] <- c(r,"upper")
    else rdat[i,] <- c(abs(r),"lower")
    
  }

  rdat[,1] <- as.double(rdat[,1])
  return(rdat)
}

wiener_likelihood <- function(x, dat) {
 if (!check_wiener_pars(x[1],x[2],x[3],x[4])) {
    return(-Inf)
  }
  ll <- vector("double", length(dat[,1]))
  for (i in 1:length(dat[,1])) {
    ll[i] <- dwiener(as.double(dat[i,1]), x[1],x[2],x[3],x[4], 
                  resp=as.character(dat[i,2]), give_log=TRUE)
  }
  return(sum(ll))
}
  
wiener_deviance <- function(x, dat) {
  -2*wiener_likelihood(x,dat)
}
wiener_bic <- function(x, dat, loss=NULL) {
  if(is.null(loss)) {
    -2*wiener_likelihood(x,dat)+4*log(length(dat[,1]))
  }
  else {
    if(is.list(dat)) {
      loss(x,dat)+length(x)*log(length(dat[[1]][,1]))
    }
    else if (is.data.frame(dat)) {
      loss(x,dat)+length(x)*log(length(dat[,1]))
    }
    else {
      stop("don't know how to handle the dat object!")
    }

  }
}
wiener_aic <- function(x, dat, loss=NULL) {
  if(is.null(loss)) {
    -2*wiener_likelihood(x,dat)+4*2 
  }
  else {
    loss(x,dat)+length(x)*2 
  }
}

# Plot function by Rainer W. Alexandrowicz
wiener_plot = function(dat)  {
  rt = as.double(dat$q)                  # response time
  rc = as.numeric(dat$resp)              # response cat: 1=up 2=lo
  dpos = try(density(rt[rc==1],from=0))  # density upper
  dneg = try(density(rt[rc==2],from=0))  # density lower
  maxt = max(pretty(max(rt)))            # overall max response time
  maxd = max(dpos$y,dneg$y)              # overall max density

  par(mar=c(0,5,0,0),mfcol=c(2,1),ask=FALSE)

  plot(dpos,xlim=c(0,maxt),ylim=c(0,maxd),las=2,lwd=2,col="green3",
      main="",ylab="",ask=FALSE)
      rug(rt[rc== 1],col="green3")
      mtext("Density of positive responses",side=2,line=4,cex=0.8)
  plot(dneg,xlim=c(0,maxt),ylim=c(maxd,0),las=2,lwd=2,col="red",
      main="",ylab="",ask=FALSE)
      mtext("Density of negative responses",side=2,line=4,cex=0.8)
      rug(rt[rc==2],col="red",side=3)
}
