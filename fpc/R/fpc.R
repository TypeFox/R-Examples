# .packageName <- "fpc"
# .First.lib <- function(lib,pkg){
#   require(cluster)
# #   if (as.numeric(R.version$major)<=1 & as.numeric(R.version$minor)<9){
# #    require(lqs)
# #    require(mva)
# #  }
#   require(MASS)
# }
#   

#
# generic functions
#
fpclusters <- function(object, ...)
  UseMethod("fpclusters")

#
# general/set distance utilities
#
#
# sseg= position of distance between i and j in smatrix/imatrix
sseg <- function (i,j){
  if (i>=j)
    out <- (i-1)*i/2+j
  else
    out <- (j-1)*j/2+i
  out
}

# Connectivity components of a graph by depth-first search
# Input: Boolean coincidence matrix, output: Vector of component numbers
con.comp <- function(comat){
  nc <- ncol(comat)
  ccn <- rep(0, times=nc) # con.comp number for each point
  fhist <- rep(FALSE, times=nc) # indicator if point had been under consideration  
  stn <- 0                  # current cc number
  pn <- 1                  # point no. to which similar objects are looked for
  while(pn>0){
    stn <- stn+1
    repeat{
      sm <- 0              # smallest new point no.
      ccn[pn] <- stn
      fhist[pn] <- TRUE
      if(nc>1)
      {
  	for(i in 2:nc)
        {
  #          cat(i, pn, ccn[i], comat[i,pn],"\n")
  	  if((ccn[i]==0) & (comat[i,pn]))
            ccn[i] <- stn
  	  if ((sm==0) & (ccn[i]==stn) & (fhist[i]==FALSE))
            sm <- i
  	} # for i
  #      cat("stn=", stn, "sm=", sm, "\n")
      } # if nc>1
      if (sm>0)
        pn <- sm
      else
        break
    } # repeat
#    print("repeat terminated")
    pn <- 0
    i <- 2
    while(i<=nc){
      if(ccn[i]==0){
        pn <- i
        i <- nc
      } # if
      i <- i+1
    } # while i
  } # while pn>0 (stn-loop)
  ccn
}    

# Intersection matrix between significant fpc groups
simmatrix <- function(fpcobj){
  sim <- c()
#  print(stn)
  for(i in 1:fpcobj$stn)
    for(j in i:fpcobj$stn){
#      cat(i,j,fpcobj$sfpc[fpcobj$sto[i]],fpcobj$sfpc[fpcobj$sto[j]])
      sim[sseg(i,j)] <- fpcobj$imatrix[sseg(fpcobj$sfpc[fpcobj$sto[i]],fpcobj$sfpc[fpcobj$sto[j]])]
    }
  sim
}

#
# fixmahal utilities
#
# Weighted ML-covariance matrix
cov.wml <- function (x, wt = rep(1/nrow(x), nrow(x)),
                     cor = FALSE, center = TRUE) 
{
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    else if (!is.matrix(x)) 
        stop("x must be a matrix or a data frame")
    if (!all(is.finite(x))) 
        stop("x must contain finite values only")
    n <- nrow(x)
    if (with.wt <- !missing(wt)) {
        if (length(wt) != n) 
            stop("length of wt must equal the number of rows in x")
        if (any(wt < 0) || (s <- sum(wt)) == 0) 
            stop("weights must be non-negative and not all zero")
        wt <- wt/s
    }
    if (is.logical(center)) {
        center <- if (center) 
            colSums(wt * x)
        else 0
    }
    else {
        if (length(center) != ncol(x)) 
            stop("length of center must equal the number of columns in x")
    }
    x <- sqrt(wt) * sweep(x, 2, center)
    cov <- (t(x) %*% x)
    y <- list(cov = cov, center = center, n.obs = n)
    if (with.wt) 
        y$wt <- wt
    if (cor) {
        sdinv <- diag(1/sqrt(diag(cov)), nrow(cov))
        y$cor <- sdinv %*% cov %*% sdinv
    }
    y
}

# inversion of cov-matrices; if singular, Eigenvalues below 1/cmax are set
# to 1/cmax.
solvecov <- function(m, cmax=1e10){
  options(show.error.messages = FALSE)
  covinv <- try(solve(m))
  if(class(covinv)!="try-error")
     coll=FALSE 
  else{ 
    p <- nrow(m)
    cove <- eigen(m, symmetric=TRUE)
    coll <- TRUE
    if (min(cove$values)<1/cmax){
      covewi <- diag(p)
      for (i in 1:p)
        if (cove$values[i]<1/cmax) covewi[i,i] <- cmax
        else covewi[i,i] <- 1/cove$values[i]
    }
    else covewi <- diag(1/cove$values,nrow=length(cove$values))
#  cat("covewi ",covewi)
    covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
  }
  options(show.error.messages = TRUE)
  out <- list(inv=covinv,coll=coll)
}

# Generation of ca by "recursion for convergence"
cmahal <- function(n,p,nmin,cmin,nc1,c1=cmin,q=1){
  cs <- rep(c1,max(n,nc1))
  for (i in (nc1-1):nmin)
    cs[i] <- cs[i+1]+q*(cs[i+1]-p)/(i+1-q)
  cs[1:(nmin-1)] <- rep(cs[nmin],nmin-1)
  if (n>nc1)
    for (i in (nc1+1):n){
      cs[i] <- cs[i-1]-q*(cs[i-1]-p)/i
      if (cs[i]<=cmin){
        cs[i] <- cmin
        i <- n
      }
    }
  cs
}

wfu <- function(md,ca,ca2,a1=1/(ca-ca2),a0=-a1*ca2){
  v <- rep(0,length(md))
  v[md<=ca] <- 1
  mind <- (md>ca) & (md<=ca2)
  v[mind] <- a1*md[mind]+a0
  v
}

# mahalanofix=vector of mahalanobis distances from all points to
# center of points indexed by gv (default=all)
mahalanofix <- function (x, n=nrow(as.matrix(x)), p=ncol(as.matrix(x)),
                         gv=rep(1, times=n), 
                         cmax=1e+10, method="ml") {
#  print("Mahalanofix")
  gv <- as.logical(gv)
  ng <- sum(gv)
  xg <- x[gv,1:p,drop=FALSE]
  if (method=="ml"){
    mg <- colMeans(xg)
    if (method=="ml")
      covg <- (ng-1)*cov(xg)/ng
    else
      covg <- cov(xg)
  }
  else{
#    if (as.numeric(R.version$major)<=1 & as.numeric(R.version$minor)<9)
#      require(lqs)
#    require(MASS)
    grob <- cov.rob(xg, method=method)
    mg <- grob$center
    covg <- grob$cov
  }
  cm <- solvecov(covg,cmax=cmax)
#  print("robmahal")
  md <- mahalanobis(x,mg,cm$inv,inverted=TRUE)
#  print("end robmahal")
  out <- list(md=md, mg=mg, covg=covg, covinv=cm$inv, coll=cm$coll)
}

# vector of mahalanobis distances from all points to
# center of points indexed by "fuzzy" gv (default=all)
mahalanofuz <- function (x, n=nrow(as.matrix(x)), p=ncol(as.matrix(x)),
                         gv=rep(1, times=n), 
                         cmax=1e+10) {
#  print("Mahalanofix")
  mg <- cov.wml(x,gv)
  covg <- mg$cov
  mg <- mg$center
  cm <- solvecov(covg, cmax=cmax)
#  print("Mahalanobis")
  md <- mahalanobis(x,mg,cm$inv,inverted=TRUE)
#  print("End Mahalanofix")
  out <- list(md=md, mg=mg, covg=covg, covinv=cm$inv, coll=cm$coll)
}

# Start configuration for fixmahal, point no. no
# startn points, cov/center of whole data
mahalconf <- function(x, no, startn, covall, plot){
  p <- ncol(x)
  n <- nrow(x)
  gv <- rep(FALSE, n)
#  print("mahalconf")
#  print(p)
#  print(n)
#  print(no)
#  print(startn)
#  print(covall)
  om <- order(mahalanobis(x,center=x[no,],cov=solvecov(covall)$inv, inverted=TRUE))
  gv[om[1:(p+1)]] <- TRUE
#  cat("conf ",gv," ",om,"\n")
  if (plot=="start" || plot=="both")
    plot(x,col=1+gv+((1:n)==no),pch=1+gv+((1:n)==no))
  if (startn>p+1)
    for (pn in (p+2):startn){
        om <- order(mahalanofix(x,n,p,gv)$md)
      j <- 0
      while (sum(gv)<pn){
        j <- j+1
        if (!gv[om[j]])
          gv[om[j]] <- TRUE
      }
    }
#   print(sum(gv))
  gv
}                

#
# fixmahal iteration
#
fpmi <- function (dat, n=nrow(as.matrix(dat)), p=ncol(as.matrix(dat)),
                  gv, ca, ca2, method="ml", plot,
                  maxit=5*n, iter=n*1e-6) {
  gva <- rep(0, times=n)
  change <- TRUE
  coll <- FALSE
  ir <- 0
  while (change) {
    ir <- ir+1
    if (method=="fuzzy")
      change <- sum((gv-gva)^2>iter)
    else
      change <- !identical(gv,gva)
    if (ir>maxit){
      change <- FALSE
      warning("Maximum number of iterations exceeded.")
    }
    gva <- gv
#  print("mahalanobis")
    if (method=="fuzzy")
      mfix <- mahalanofuz(x=dat, n=n, p=p, gv=gv)
    else
      mfix <- mahalanofix(x=dat, n=n, p=p, gv=gv, method=method)
#  print(ca)
    md <- mfix$md
    mg <- mfix$mg
    if (mfix$coll) coll=TRUE
    covg <- mfix$covg
#  print(md)
    if (method=="fuzzy"){
      gv <- wfu(md, ca, ca2)
#      print(gv)
      cn <- ca
#     plot(x,col=1+gv)
    }
    else{
      ng <- sum(gv)
      if (length(ca)==1)
        cn <- ca
      else
        cn <- ca[ng]
#    print(ng)
#    cat("Target function = ",
#      det(covg*(ng-1)/ng)*prod(exp((2-ca)/(1:ng))), "\n")
      gv <- md<=cn
    }
#  print(gv)
     if (plot=="iteration" || plot=="both"){
       if (method=="fuzzy")
         plotgv <- as.integer(gv>0.99)+
           3*as.integer(gv>0.01 & gv<=0.99)
       else
         plotgv <- gv
       plot(dat,col=1+plotgv,pch=1+plotgv)
     }
  } # while change
  out <- list(mg=mg, covg=covg, md=md, gv=gv, coll=coll, method=method, ca=cn)
  out
}

#
# fixmahal main and output
#
fixmahal <- function (dat, n=nrow(as.matrix(dat)), p=ncol(as.matrix(dat)), 
                      method="fuzzy", cgen="fixed", 
                      ca=NA, ca2=NA,
                      calpha=ifelse(method=="fuzzy",0.95,0.99),
                      calpha2=0.995,
                      pointit=TRUE, subset=n,
                      nc1=100+20*p,
                      startn = 18+p, mnc = floor(startn/2), 
                      mer=ifelse(pointit,0.1,0), distcut=0.85,  
                      maxit=5*n, iter=n*1e-5, 
                      init.group=list(), 
                      ind.storage=TRUE, countmode=100, 
                      plot="none"){

# ca2, calpha2 are ignored unless method="fuzzy"
# Initializations, WDC

# print("fixmahal")
  if (startn>n) startn <- n
  dat <- as.matrix(dat)
  if (method=="fuzzy"){
    if (is.na(ca)){
      ca <- qchisq(calpha,p)
      ca2 <- qchisq(calpha2,p)
    }
  }
  else{
    if (cgen=="fixed"){
      if (is.na(ca))
        ca <- qchisq(calpha,p)
      else
        calpha <- pchisq(ca,p)
    }
    else{
      if (is.na(ca))
        ca <- qchisq(calpha,p)
      ca <- cmahal(n,p,nmin=startn,cmin=ca,nc1=nc1,q=1)
    }
    ca2 <- 0
  }  
  dat <- dat[1:n,1:p,drop=FALSE]
  tsc <- 0          # number of too small FPCs
  ncoll <- 0        # number of iterations leading to collinear data
  gv <- rep (1, times=n)
# print("fpmi")
  fpc <- fpmi(dat, n, p, gv, ca, ca2, method, plot, maxit, iter)   # WDC-iteration
# print("1Ende")
  if (cgen=="auto")
    cv <- c(fpc$ca)
  imatrix <- c(sum(fpc$g)) # matrix of FPC size and intersections
  smatrix <- c(sum(fpc$g)) # matrix of FPC size and similarities
  nc <- 1                  # number of FPCs
  clist <- list(fpc$mg)  # list of FPC locations
  vlist <- list(fpc$covg)   # list of FPC covariance matrices
  nfound <- c(1)           # times found per FPC
  expectratio <- c()       # nfound/size
  if (ind.storage)
    glist <- list(fpc$g)   # list of FPC indicator vectors

# Standard iterations

  if (pointit){
    if(subset < n)
      itpoints <- sample(1:n,subset)
    else{
      itpoints <- 1:n
      subset <- n
    }
    if (method=="fuzzy" || method=="classical" || method=="ml")
      covall <- cov(dat)
    else
      covall <- cov.rob(dat, method=method)$cov
    for (ni in 1:subset){
      i <- itpoints[ni]
      if(countmode*round(ni/countmode)==i)
        cat("Iteration run no. ", ni, " of ",subset, "\n")
# print(covall)
      gv <- mahalconf(dat,i,startn,covall=covall, plot)
#  cat("i= ",i," gv=  ",(1:366)[gv]," \n")
# print("iteration")
  #      cat(sum(gv[1:50]),sum(gv[51:100]))
      fpc <- fpmi(dat, n, p, gv, ca, ca2, method, plot, maxit, iter)
#   print(fpc$mg)
  # print("neu?")
      neu <- TRUE
      j <- 1
      while(j<=nc) {
        if (sum((clist[[j]]-fpc$mg)^2)<1e-6 && sum((vlist[[j]]-fpc$covg)^2)<1e-5){
          neu <- FALSE
          cnum <- j
        }          # if j found
        j <- j+1
      }            # while j
      ng <- sum(fpc$g)
      if (neu & (ng>=mnc)){
        nc <- nc +1
        nfound[nc] <- 1
        clist[[nc]] <- fpc$mg
        vlist[[nc]] <- fpc$covg
        if (cgen=="auto")
          cv <- c(cv,fpc$ca)
# cat("ind.storage= ",ind.storage,  "nc= ", nc, "\n")
        if (ind.storage){
          glist[[nc]] <- fpc$g
          j <- 1
          while( j<=(nc-1)) {
            imatrix[sseg(nc,j)] <- sum(pmin(fpc$g,glist[[j]]))
# cat("   (storage T) smatrix: ",nc,j,smatrix[sseg(nc,j)],"\n")
            j<- j+1
          }        # for j
        }          # if ind.storage
        else{
          for(j in 1:(nc-1)){
            if (method=="fuzzy"){
              mah <- mahalanobis(dat,center=clist[[j]],
                           cov=solvecov(vlist[[j]])$inv,inverted=TRUE)
              glistj <- wfu(mah, ca, ca2)
            }
            else{
              if (cgen=="fixed"){
                  glistj <- (mahalanobis(dat,center=clist[[j]],
                           cov=solvecov(vlist[[j]])$inv,inverted=TRUE) <= ca)
              }
              else{
                  glistj <- (mahalanobis(dat,center=clist[[j]],
                           cov=solvecov(vlist[[j]])$inv,inverted=TRUE) <= cv[j])
              }
            }
            imatrix[sseg(nc,j)] <- sum(pmin(fpc$g,glistj))
# cat("   smatrix: ",nc,j,smatrix[sseg(nc,j)],"\n")
          }       # for j
        }          # else (!ind.storage)
        imatrix[sseg(nc,nc)] <- sum(fpc$g)
#        cat("Neu: Punkt ",i,", Cluster mit ",sum(fpc$g)," Punkten. \n")
      }            # if neu & ng>=mnc
      if (ng<mnc)
        tsc <- tsc+1
      if (!neu){
        nfound[cnum] <- nfound[cnum]+1
      }
      if (neu & fpc$coll)
        ncoll <- ncoll+1
    }              # for i
  }                # if pointit
# Iterations with init.group

#  print(length(init.group))
  grfpc <- FALSE
  if(length(init.group)>0){ 
    i <- 1
    grfpc <- rep(0, times=length(init.group))
    while(i<=length(init.group)){
      gv <- init.group[[i]]  
#      print(sum(gv))
#      print(length(gv))
      fpc <- fpmi(dat, n, p, gv, ca, ca2, method, plot, maxit, iter)
#      print("fpc")
      neu <- TRUE
      j <- 1
      while(j<=nc) {
        if (sum((clist[[j]]-fpc$mg)^2)<1e-6 && sum((vlist[[j]]-fpc$covg)^2)<1e-5){
          neu <- FALSE
          cnum <- j
        }          # if j found
        j <- j+1
      }            # while j
      ng <- sum(fpc$g)
      if (neu & (ng>=mnc)){
# print("neu! smatrix")
        nc <- nc +1
        nfound[nc] <- 1
        clist[[nc]] <- fpc$mg
        vlist[[nc]] <- fpc$covg
        if (cgen=="auto")
          cv <- c(cv,fpc$ca)
# cat("ind.storage= ",ind.storage,  "nc= ", nc, "\n")
        if (ind.storage){
          glist[[nc]] <- fpc$g
          j <- 1
          while( j<=(nc-1)) {
            imatrix[sseg(nc,j)] <- sum(pmin(fpc$g,glist[[j]]))
# cat("   (storage T) smatrix: ",nc,j,smatrix[sseg(nc,j)],"\n")
            j<- j+1
          }        # for j
        }          # if ind.storage
        else{
          for(j in 1:(nc-1)){
            if (method=="fuzzy"){
              mah <- mahalanobis(dat,center=clist[[j]],
                           cov=solvecov(vlist[[j]])$inv,inverted=TRUE)
              glistj <- wfu(mah, ca, ca2)
            }
            else{
              if (cgen=="fixed"){
                  glistj <- (mahalanobis(dat,center=clist[[j]],
                           cov=solvecov(vlist[[j]])$inv,inverted=TRUE) <= ca)
              }
              else{
                  glistj <- (mahalanobis(dat,center=clist[[j]],
                           cov=solvecov(vlist[[j]])$inv,inverted=TRUE) <= cv[j])
              }
            }
            imatrix[sseg(nc,j)] <- sum(pmin(fpc$g,glistj))
# cat("   smatrix: ",nc,j,smatrix[sseg(nc,j)],"\n")
          }       # for j
        }          # else (!ind.storage)
        imatrix[sseg(nc,nc)] <- sum(fpc$g)
      }            # if neu & ng>=mnc
      if (ng<mnc)
        tsc <- tsc+1
      if (!neu)
        nfound[cnum] <- nfound[cnum]+1
      else
        cnum <- nc
      if (neu & fpc$coll)
        ncoll <- ncoll+1
      grfpc[i] <- cnum
      i <- i+1
    }             # while i
  }               # if init.group

# Find structures

# print(nc)
 if(nc>1){
    for (i in 1:(nc-1)){
      smatrix[sseg(i,i)] <- imatrix[sseg(i,i)]
      for (j in (i+1):nc){
        smatrix[sseg(i,j)] <- 2*imatrix[sseg(i,j)]/(imatrix[sseg(i,i)]+
                              imatrix[sseg(j,j)])
#        cat("smatrix: ",i,j,smatrix[sseg(i,j)],"\n")
      } # for j
    }
    smatrix[sseg(nc,nc)] <- imatrix[sseg(nc,nc)]
    comat <- matrix(nrow=nc,ncol=nc)
    for (i in 1:(nc-1))
      for(j in (i+1):nc)
        comat[i,j] <- comat[j,i] <- (smatrix[sseg(i,j)]>=distcut)
  } #if nc>1
  else
    comat <- matrix(1,1)
#  print(comat)
  struc <- con.comp(comat)
#  print(struc)
  stn <- max(struc)
  rm(comat)

# print("description of structures")

  tf <- rep(0, times=stn)   # structure: times found
  maxf <- rep(0, times=stn) # structure, rep. fpc: expectratio
  sfpc <- rep(0, times=stn) # structure: representative fpc
  ser <- rep(0, times=stn) # structure: expectation ratio
  for(i in 1:nc){
    expectratio[i] <- nfound[i] / imatrix[sseg(i,i)]
    tf[struc[i]] <- tf[struc[i]] + nfound[i]
    if(expectratio[i]==maxf[struc[i]])
      if(imatrix[sseg(i,i)]<imatrix[sseg(sfpc[struc[i]],sfpc[struc[i]])])
        sfpc[struc[i]] <- i
    if(expectratio[i]>maxf[struc[i]]){
        sfpc[struc[i]] <- i
        maxf[struc[i]] <- expectratio[i]
    }
  } # for i
  for (i in 1:stn)
    ser[i] <- (tf[i]*n)/(imatrix[sseg(sfpc[i],sfpc[i])]*subset)
  skc <- stn
  stn <- sum(ser>=mer)
  skc <- skc-stn
#  print(stn)
  
# Output
  cvec <- ca
  if (cgen=="auto")
    ca <- cv
  if (!ind.storage)
    glist=FALSE
  out <- list(nc=nc, g=glist, means=clist, covs=vlist,
                nfound=nfound, er=expectratio,
                tsc=tsc, ncoll=ncoll, skc=skc, grto=grfpc,
                imatrix=imatrix, smatrix=smatrix,
                stn=stn, stfound=tf, ser=ser, 
                sfpc=sfpc, ssig=sfpc[ser>=mer], sto=order(-ser), 
                struc=struc, 
                n=n, p=p, method=method, cgen=cgen, 
                ca=ca, ca2=ca2, cvec=cvec, calpha=calpha,
                pointit=pointit, subset=subset, 
                mnc=mnc, startn=startn, mer=mer, distcut=distcut)
  class(out) <- "mfpc"
  out   
}

summary.mfpc <- function(object, ...){
#  print("Beginn reducemahal")
  clist <- list(0)
  vlist <- list(0)
  expectratio <- c()
  tf <- c()
  sn <- c()
  ca <- object$ca
  method <- object$method
  strn <- object$stn
  if (strn>0)
    for(i in 1:strn){
#    cat(i, object$coefs[[object$sfpc[object$sto[i]]]], "\n")
      clist[[i]] <- object$means[[object$sfpc[object$sto[i]]]]
#    cat(i, object$vars[[object$sfpc[object$sto[i]]]], "\n")
      vlist[[i]] <- object$covs[[object$sfpc[object$sto[i]]]]
      tf[i] <- object$stfound[object$sto[i]]
      sn[i] <- object$imatrix[sseg(object$sfpc[object$sto[i]],object$sfpc[object$sto[i]])]
      expectratio[i] <- object$ser[object$sto[i]]
      if (object$cgen=="auto")
        ca[i] <- object$ca[object$sfpc[object$sto[i]]]
    }
  tskip <- sum(object$nfound)-sum(tf)
  if (strn>0)
    sim <- simmatrix(object)
  else
    sim <- NULL
  out <- list(means=clist, covs=vlist, stn=strn, stfound=tf, sn=sn, 
              ser=expectratio, 
              tskip=tskip, skc=object$skc, tsc=object$tsc, sim=sim, 
              ca=ca, ca2=object$ca2, calpha=object$calpha,
              mer=object$mer, mnc=object$mnc,
              method=method, cgen=object$cgen,
              pointit=object$pointit)
  class(out) <- "summary.mfpc"
  out
}

# clusters: List of indicator vectors of representative 
# fpcs of significant structures
# in order of stfound.
# Specify x if ind.storage in fixreg was F
fpclusters.mfpc <- function(object, dat=NA, ca=object$ca, p=object$p, ...){
  glist <- list()
  if (object$method=="fuzzy"){
    ca <- object$ca
    ca2 <- object$ca2
  }
# cat("stn= ",object$stn,"\n")
  if (object$stn>0)
    for(i in 1:object$stn){
      if (identical(object$g,FALSE)){
        mc <- object$means[[ object$sfpc[object$sto[i]] ]]
        if (p>1){
          cm <- solvecov(object$covs[[ object$sfpc[object$sto[i]] ]])
          if (object$method=="fuzzy")
            glist[[i]] <- wfu(mahalanobis(dat,mc,cm$inv,inverted=TRUE), ca, ca2)
          if (object$cgen=="fixed")
            glist[[i]] <- (mahalanobis(dat,mc,cm$inv,inverted=TRUE)<=ca)
          else
            glist[[i]] <- (mahalanobis(dat,mc,cm$inv,inverted=TRUE)<=ca[object$sfpc[object$sto[i]]] )
        }
        else{
          cm <- object$covs[[ object$sfpc[object$sto[i]] ]]
          if (object$cgen=="fixed")
            glist[[i]] <- as.vector((dat-mc)^2)<=as.double(ca*cm)
          else
            glist[[i]] <- as.vector((dat-mc)^2)<=as.double(ca[object$sfpc[object$sto[i]]]*cm)
        }
      }   # if g==F
      else
        glist[[i]] <- object$g[[object$sfpc[object$sto[i]]]]
  # print(i)
    } # for i
  else
    warning("No FPC was found often enough!")
  glist
}

# Bhattacharyya coordinates of FPC no. no
# bw: black/white
plot.mfpc <- function(x, dat, no, bw=FALSE,
                      main=c("Representative FPC No. ",no),
                      xlab=NULL, ylab=NULL,
                      pch=NULL, col=NULL, ...){
  n <- x$n
  p <- x$p
  ca <- x$ca
  sumobj <- summary(x)
  if (x$cgen=="auto")
    ca <- sumobj$ca[no]
  mg <- sumobj$means[[no]]
  covg <- sumobj$covs[[no]]
#  print("clusters")
  if (x$method=="fuzzy"){
    g <- fpclusters(x,dat)[[no]] 
    gv <- as.integer(g>0.5)
    gvp <- as.integer(g>0.99)+
      3*as.integer(g>0.01 & g<=0.99)
  }
  else
    gvp <- gv <- as.integer(fpclusters(x,dat)[[no]])
  if (n==sum(gv)){
    cat("Cluster equals whole dataset => No cluster plot is done.\n")
  }
  else{
    if (is.null(pch))
      pch <- if (bw) gvp+1 else 1
    if (is.null(col))
      col <- if (bw) 1 else gvp+1
#  print("plot")
    if (p==1){
        if (is.null(xlab))
          xlab <- "Index"
        if (is.null(ylab))
          ylab <- deparse(substitute(dat))
        plot(dat,col=col, pch=pch, main=main, xlab=xlab, ylab=ylab, ...)
        abline(c(mg,0))
        abline(c(mg-sqrt(ca*covg),0), lty="dotted")       
        abline(c(mg+sqrt(ca*covg),0), lty="dotted")
    }
    if (p==2){
        xy <- xy.coords(dat,NULL)
        if (is.null(xlab))
          xlab <- xy$xlab
        if (is.null(ylab))
          ylab <- xy$ylab
        plot(dat,col=col, pch=pch, main=main, xlab=xlab, ylab=ylab, ...)
        points(as.vector(mg[1]),as.vector(mg[2]),pch="M",
                 col=ifelse(bw, col, "blue"), ...)
        cge <- eigen(covg, symmetric=TRUE)
        poly <- rep(0, times=200)
        dim(poly) <- c(100,2)
        for(i in 1:100)
          poly[i,] <- (mg + c(sin(i*2*pi/100), cos(i*2*pi/100)) %*% 
            (t(cge$vectors) * sqrt(cge$values))* sqrt(ca))
        polygon(poly[,1], poly[,2], border="black")  
    }
    if (p>2){
      gc1 <- batcoord(dat,gv)
        if (is.null(xlab))
          xlab <- "Discriminant projection 1"
        if (is.null(ylab))
          ylab <- "Discriminant projection 2"
        plot(gc1$proj[,1:2],col=col, pch=pch,
             main=main, xlab=xlab, ylab=ylab, ...)
        points(as.vector((mg %*% gc1$units)[1]),
               as.vector((mg %*% gc1$units)[2]), pch="M",
                 col=ifelse(bw, col, "blue"), ...)
    }
  }
  invisible()
}

print.mfpc <- function(x, ...){
  cat("Mahalanobis Fixed Point Cluster object\n")  
  cat(x$stn," representative stable fixed point clusters\n") 
  cat(" of totally ",x$nc," found fixed point clusters.\n")
  invisible(x)
}

# Summary output of mfpc objects; choose fpcobj <- summary(fpcobj) to
# replace fpcobj by its reduction    
print.summary.mfpc <- function(x, maxnc=30, ...){
#  print("Beginn Summary")
#  print(class(x))
  minnc <- min(x$stn,maxnc)
  cat("  *  Mahalanobis Fixed Point Clusters  *\n\n")
  cat("Often a clear cluster in the data leads to several similar FPCs.\n")
  cat("The summary shows the representative FPCs of groups of similar FPCs.\n\n")
  cat("Method ",x$method," was used.\n")
  cat("Number of representative FPCs: ", x$stn, "\n\n")
  cat("FPCs with less than ",x$mnc," points were skipped.\n")
  if (x$mer>0) 
    cat("FPCs with ratio of times found to number of points less than ",
        x$mer," were skipped.\n")
  cat(x$tsc+x$tskip, " iteration runs led to ",x$skc+x$tsc," skipped clusters.\n")
  if (x$stn>maxnc)
    cat("Warning! Only ",maxnc," clusters are displayed. \nSpecify maxnc in call of summary.mfpc if you want enlarged display or \nrun fixmahal with larger ca, calpha, mnc or mer to get fewer clusters. \n")
  if (x$method=="fuzzy")
    cat("  Weight 1 for r^2<= ",x$ca," weight 0 for r^2> ",
    x$ca2," \n")   
  if (x$cgen=="fixed")
    cat("  Constant ca= ", x$ca,
      " corresponding to alpha= ",x$calpha,"\n")
  cat("\n")
  if (x$stn==0)
    cat("No FPC was found often enough.\n")
  else{
    for(i in 1:minnc){
      cat(" FPC ",i, "\n")
      cat("  Times found (group members): ",x$stfound[i], "\n")
      if (x$pointit)
        cat("  Ratio to size: ",x$ser[i], "\n")
      cat("  Mean:\n")
      print(x$means[[i]])
      cat("  Covariance matrix:\n")
      print(x$covs[[i]])
      if (x$cgen=="auto" && x$method!="fuzzy")
        cat("  Constant ca= ",x$ca[[i]],"\n")
      cat("  Number of points (sum of weights): ",x$sn[[i]], "\n\n")
    }
    cat("Number of points (rounded weights) in intersection of representative FPCs\n")
    sm <- rep(0,times=minnc^2)
    dim(sm) <- c(minnc,minnc)
    for(i in 1:minnc)
      for(j in 1:minnc)
        sm[i,j] <- round(x$sim[sseg(i,j)])
    print(sm)
  }
  invisible(x)
}











