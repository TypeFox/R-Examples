# rewrite crr in 'cmprsk' package
# all the same, except no censoring group option
# all the C functions called are in crrc.c
# completed: 12/21/2008
# author: Bingqing Zhou

crrc <- function(ftime,fstatus,cov1,cov2,tf,cluster, cengroup,failcode=1,
          cencode=0, subset,na.action=na.omit,gtol=1e-6,maxiter=10,init) {

  d <- data.frame(ftime=ftime,fstatus=fstatus,
        cengroup=if (missing(cengroup)) rep(1,length(fstatus)) else cengroup,
        cluster=if (missing(cluster)) 1:length(fstatus) else cluster)

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
  if (nrow(d) != tmp) cat(format(tmp-nrow(d)),
        'cases omitted due to missing values\n')

  d = d[order(d$ftime),]
  ftime = d$ftime
  cenind = ifelse(d$fstatus==cencode,1,0)
  fstatus = ifelse(d$fstatus==failcode,1,2*(1-cenind))
  ucl = sort(unique.default(d$cluster))
  cluster = match(d$cluster,ucl)
  ncl = length(ucl)  
  ucg = sort(unique.default(d$cengroup))
  cengroup = match(d$cengroup,ucg)
  ncg = length(ucg)

### want censring dist km at ftime-
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
    cov1 = as.matrix(d[,(1:np1)+4])
    cov2 = 0
    tfs = matrix(0,length(ft),1)
  } else if (np1 == 0) {
    cov2 = as.matrix(d[,(1:np2)+4+np1])
    cov1 = 0
    tfs = as.matrix(tf(ft))
  } else {
    cov1 = as.matrix(d[,(1:np1)+4])
    cov2 = as.matrix(d[,(1:np2)+4+np1])
    tfs = as.matrix(tf(ft))
  }

### start of nr
  if (missing(init)) b = rep(0,np)
  else b = init
  stepf = .5

  for (ll in 0:maxiter) {
    # only for 1 cengroup
    z <- .C("crrfsvoo", as.double(ftime),as.integer(fstatus),
           as.integer(length(ftime)), as.double(cov1), as.integer(np1),
           as.double(cov2),as.integer(np2), as.double(tfs), as.integer(nf),
           as.double(uuu), as.double(b), double(1), double(np), 
           double(np*np),PACKAGE="crrSC")[12:14]

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
    fbn <- - .C("crrfoo", as.double(ftime), as.integer(fstatus),
              as.integer(length(ftime)), as.double(cov1), as.integer(np1),
              as.double(cov2), as.integer(np2), as.double(tfs),
              as.integer(nf), as.double(uuu), as.double(bn), double(1),PACKAGE="crrSC")[[12]]

# backtracking loop
    i = 0
    while (is.na(fbn) || fbn> -z[[1]]-(1e-4)*sum(sc*z[[2]])) {
      i = i+1
      sc = sc*stepf
      bn = b+sc
    fbn <- - .C("crrfoo", as.double(ftime), as.integer(fstatus),
                as.integer(length(ftime)), as.double(cov1), as.integer(np1),
                as.double(cov2), as.integer(np2), as.double(tfs),
                as.integer(nf), as.double(uuu), as.double(bn), double(1),PACKAGE="crrSC")[[12]]
      if (i>20) break
    }
    
    if (i>20) {
      converge <- FALSE
      break
    }
    b = c(bn)
  }

   vout <- .C("crrvvoo", as.double(ftime),as.integer(fstatus),
          as.integer(length(ftime)), as.double(cov1), as.integer(np1),
          as.double(cov2), as.integer(np2),as.double(tfs),
          as.integer(nf), as.integer(cluster),as.integer(ncl), as.double(uuu),  
          as.double(b), double(np*np), double(np*np),PACKAGE="crrSC")[14:15]
          
  dim(vout[[1]]) = dim(vout[[2]]) = c(np,np)
    
  h = solve(vout[[1]])
  v = h %*% vout[[2]] %*% t(h)

  z = list(coef=b,loglik=z[[1]],score=z[[2]],inf=matrix(z[[3]], np, np), var=v,
      ftime=ft, tfs=as.matrix(tfs), converged=converge)
  class(z) = 'crr'
  z
}
