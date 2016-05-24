# description: goodness of fit test for proportional subdistribution hazards (psh)
# functions: "pshtest" is the score test for gof of psh
#            "crr.test" performs the key calculations in pshtest
#            "mycrr" is a rewritten of the "crr" function in "cmprsk' package (no censoring group)
# Some examples in the end of the file 
# all the C functions called are in crrc.c (completed: 12/21/2008, updated June 2013)
# reference: 'cmprsk' package by Robert Gray
# updated: 6/19/2013
# author: Bingqing Zhou

#library(survival)

mycrr <- function(ftime,fstatus,cov1,cov2,tf,cengroup,failcode=1,cencode=0,
         subset,na.action=na.omit,gtol=1e-6,maxiter=10,init) {

  d <- data.frame(ftime=ftime,fstatus=fstatus,
        cengroup=if (missing(cengroup)) rep(1,length(fstatus)) else cengroup)

  if (!missing(cov1)) {
    cov1 = as.matrix(cov1)
    np1 = ncol(cov1)
    d = cbind(d,cov1)
  } else {np1 = 0}

  if (!missing(cov2)) {
    cov2 = as.matrix(cov2)
    np2 = ncol(cov2)
    d = cbind(d,cov2)
  } else {np2 = 0}
  
   np = np1 + np2

  if (!missing(subset)) d = d[subset,]
  tmp = nrow(d)
  d = na.action(d)
  if (nrow(d) != tmp) cat(format(tmp-nrow(d)),'cases omitted due to missing values\n')

  d = d[order(d$ftime),]
  ftime = d$ftime
  cenind = ifelse(d$fstatus==cencode,1,0)
  fstatus = ifelse(d$fstatus==failcode,1,2*(1-cenind))
  ucg = sort(unique.default(d$cengroup))
  cengroup = match(d$cengroup,ucg)
  ncg = length(ucg)

### want censoring dist km at ftime-
  uuu = matrix(0,nrow=ncg,ncol=length(ftime))
  for (k in 1:ncg) {
    u = do.call('survfit',list(formula=Surv(ftime,cenind)~1,
            data=data.frame(ftime,cenind,cengroup),subset=cengroup==k))
    u = summary(u,times=sort(ftime*(1-.Machine$double.eps)))
    uuu[k,1:length(u$surv)] = u$surv
  }

  ft = sort(ftime[fstatus==1])
  nf = length(ft)
  if (np2 == 0) {
    cov1 = as.matrix(d[,(1:np1)+3])
    cov2 = 0
    tfs = matrix(0,length(ft),1)
  } else if (np1 == 0) {
    cov2 = as.matrix(d[,(1:np2)+3+np1])
    cov1 = 0
    tfs = as.matrix(tf(ft))
  } else {
    cov1 = as.matrix(d[,(1:np1)+3])
    cov2 = as.matrix(d[,(1:np2)+3+np1])
    tfs = as.matrix(tf(ft))
  }

### start of nr
  if (missing(init)) b = rep(0,np)
  else b = init
  stepf = .5

  for (ll in 0:maxiter) {
    # only for 1 cengroup
    z <- .C("crrfsvo", as.double(ftime),as.integer(fstatus),
           as.integer(length(ftime)), as.double(cov1), as.integer(np1),
           as.double(cov2),as.integer(np2), as.double(tfs), as.integer(nf),
           as.double(uuu), as.double(b), double(1), double(np), double(np*np))[12:14]

    if (max(abs(z[[2]])*pmax(abs(b),1)) < max(abs(z[[1]]),1)*gtol) {
      converge = TRUE
      break
    }

    if (ll==maxiter) {
      converge = FALSE
      break
    }

    h = z[[3]]
    dim(h) = c(np,np)
    sc = solve(h,z[[2]])
    bn = b + sc
    fbn <- - .C("crrfo", as.double(ftime), as.integer(fstatus),
              as.integer(length(ftime)), as.double(cov1), as.integer(np1),
              as.double(cov2), as.integer(np2), as.double(tfs),
              as.integer(nf), as.double(uuu), as.double(bn), double(1))[[12]]

# backtracking loop
    i = 0
    while (is.na(fbn) || fbn> -z[[1]]-(1e-4)*sum(sc*z[[2]])) {
      i = i+1
      sc = sc*stepf
      bn = b+sc
    fbn <- - .C("crrfo", as.double(ftime), as.integer(fstatus),
                as.integer(length(ftime)), as.double(cov1), as.integer(np1),
                as.double(cov2), as.integer(np2), as.double(tfs),
                as.integer(nf), as.double(uuu), as.double(bn), double(1))[[12]]
      if (i>20) break
    }
    
    if (i>20) {
      converge <- FALSE
      break
    }
    b = c(bn)
  }

   vout <- .C("crrvvo", as.double(ftime),as.integer(fstatus),
          as.integer(length(ftime)), as.double(cov1), as.integer(np1),
          as.double(cov2), as.integer(np2),as.double(tfs),
          as.integer(nf), as.double(uuu),  as.double(b),
          double(np*np), double(np*np))[12:13]
          
  dim(vout[[1]]) = dim(vout[[2]]) = c(np,np)
    
  h = solve(vout[[1]])
  v = h %*% vout[[2]] %*% t(h)

  z = list(coef=b,loglik=z[[1]],score=z[[2]],inf=matrix(z[[3]], np, np), var=v,
      ftime=ft, tfs=as.matrix(tfs), converged=converge)
  class(z) = 'crr'
  z
}


#-------------------------

mycrr.test <- function(ftime,fstatus,cov1,cov2,tf,cengroup,failcode=1,cencode=0,
         subset,na.action=na.omit,gtol=1e-6,maxiter=10,init) {

  d <- data.frame(ftime=ftime,fstatus=fstatus,
        cengroup=if (missing(cengroup)) rep(1,length(fstatus)) else cengroup)

  if (!missing(cov1)) {
    cov1 = as.matrix(cov1)
    np1 = ncol(cov1)
    d = cbind(d,cov1)
  } else {np1 = 0}

  if (!missing(cov2)) {
    cov2 = as.matrix(cov2)
    np2 = ncol(cov2)
    d = cbind(d,cov2)
  } else {np2 = 0}
  
   np = np1 + np2

  if (!missing(subset)) d = d[subset,]
  tmp = nrow(d)
  d = na.action(d)
  if (nrow(d) != tmp) cat(format(tmp-nrow(d)),'cases omitted due to missing values\n')

  d = d[order(d$ftime),]
  ftime = d$ftime
  cenind = ifelse(d$fstatus==cencode,1,0)
  fstatus = ifelse(d$fstatus==failcode,1,2*(1-cenind))
  ucg = sort(unique.default(d$cengroup))
  cengroup = match(d$cengroup,ucg)
  ncg = length(ucg)

### want censoring dist km at ftime-
  uuu = matrix(0,nrow=ncg,ncol=length(ftime))
  for (k in 1:ncg) {
    u = do.call('survfit',list(formula=Surv(ftime,cenind)~1,
            data=data.frame(ftime,cenind,cengroup),subset=cengroup==k))
    u = summary(u,times=sort(ftime*(1-.Machine$double.eps)))
    uuu[k,1:length(u$surv)] = u$surv
  }

  ft = sort(ftime[fstatus==1])
  nf = length(ft)
  if (np2 == 0) {
    cov1 = as.matrix(d[,(1:np1)+3])
    cov2 = 0
    tfs = matrix(0,length(ft),1)
  } else if (np1 == 0) {
    cov2 = as.matrix(d[,(1:np2)+3+np1])
    cov1 = 0
    tfs = as.matrix(tf(ft))
  } else {
    cov1 = as.matrix(d[,(1:np1)+3])
    cov2 = as.matrix(d[,(1:np2)+3+np1])
    tfs = as.matrix(tf(ft))
  }

  b = init
    z <- .C("crrfsvo", as.double(ftime),as.integer(fstatus),
           as.integer(length(ftime)), as.double(cov1), as.integer(np1),
           as.double(cov2),as.integer(np2), as.double(tfs), as.integer(nf),
           as.double(uuu), as.double(b), double(1), double(np), double(np*np))[12:14]

   v <- .C("crrvvo", as.double(ftime),as.integer(fstatus),
          as.integer(length(ftime)), as.double(cov1), as.integer(np1),
          as.double(cov2), as.integer(np2),as.double(tfs),
          as.integer(nf), as.double(uuu),  as.double(b),
          double(np*np), double(np*np))[12:13]
          
    dim(v[[2]]) <- dim(v[[1]]) <- c(np, np)

    out <- list(coef = b, loglik = -z[[1]], score = -z[[2]], inf = matrix(z[[3]], 
        np, np), var = v[[2]], tfs = as.matrix(tfs))
    class(out) <- "crr"
    out
}

#------------
psh.test <- function(time, fstatus, z, D=c(1,1), tf=function(x) cbind(x,x^2), 
  init){
 # time: failure time, n*1
 # fstatus: failure status =0 if censored, n*1
 # z: covariates n*np
 # D: components of z that are tested for time-varying effect, nb*1
 # tf: functions of t for z being tested on the same location, nb
      z <- as.matrix(z)
      np1 = ncol(z)
      np2 = ncol(as.matrix(D))
      np = np1 + np2
    
      weout <- mycrr(ftime=time, fstatus=fstatus, cov1=z)
      weout2 <- mycrr.test(ftime=time,fstatus=fstatus,cov1=z,cov2=z %*% D,
           tf=tf, init= c(weout$coef,rep(0,np2)))

      v = matrix(weout2$inf[1:np1,1:np1], ncol=np1)
      vg = matrix(weout2$inf[(np1+1):np,1:np1],ncol=np1)
      vgg = matrix(weout2$inf[(np1+1):np,(np1+1):np],ncol=np2)
      v2 = matrix(weout2$var[1:np1,1:np1], ncol=np1)
      v2g = matrix(weout2$var[(np1+1):np,1:np1], ncol=np1)
      v2gg = matrix(weout2$var[(np1+1):np,(np1+1):np],ncol=np2)
      A = v2gg + vg %*% solve(v) %*% v2 %*% solve(v) %*% t(vg) -
          2 * vg %*% solve(v) %*% t(v2g)
      utheta = weout2$score[(np1+1):np]
      T = t(utheta) %*% solve(A) %*% utheta
      p = 1- pchisq(T, df= np2)

      out <- cbind(mean(fstatus == 0),mean(fstatus == 1), 
        round(T,4), np2, round(p,4))
      dimnames(out)[[2]] <- c('% cen','% cause 1', 'Test Statistic', 'd.f.','p-value')
      return(out)
}


# end of the file

