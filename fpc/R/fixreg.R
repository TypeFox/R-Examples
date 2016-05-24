#
# fixreg utilities
#
# Generation of ca by formula
can <- function (n,p){
    ca <- 3+33/(n*2^(-(p-1)/2))^(1/3)+2900000/(n*2^(-(p-1)/2))^3
  ca
}  
      
# clusexpect= Expectation of times found for a cluster of given size under
# assumption that it will be always found with only its points and never else.
clusexpect <- function (n, p, cn, ir) {
  prob <- choose (cn, p+2)/choose (n, p+2)
  result <- ir*prob
  result
}
 
# Calculates iteration number such that theoretical probability that cluster
# with size cn is found at least mtf times exceeds prob.   
itnumber <- function (n, p, cn, mtf, prob=0.95, maxir=20000){
  minir <- 1
  repeat
  {
     if (qbinom(1-prob,minir,clusexpect(n,p,cn,1))>=mtf){
       ir <- minir
       break
     }
     if (qbinom(1-prob,maxir,clusexpect(n,p,cn,1))<mtf){
       ir <- maxir
       warning("Iteration number too large. maxir chosen.")
       break
     }
     ir <- round((minir+maxir)/2)
     if (qbinom(1-prob,ir,clusexpect(n,p,cn,1))<mtf)
       minir <- ir
     else
       maxir <- ir
     if (maxir <= minir+1)
       break
  }     # repeat
  ir
}

# Minimum size of FPCs that will be found at least mtf times with probability
# prob under ir iteration runs.
minsize <- function(n, p, ir, mtf, prob=0.5){
  mincn <- p+2
  maxcn <- n-1
  repeat
  {
     if (qbinom(1-prob,ir,clusexpect(n,p,mincn,1))>=mtf){
       cn <- mincn
       break
     }
     if (qbinom(1-prob,ir,clusexpect(n,p,maxcn,1))<mtf){
       cn <- maxcn
       warning("Too few iteration runs. Informative FPCs are not probable.")
       break
     }
     cn <- round((mincn+maxcn)/2)
     if (qbinom(1-prob,ir,clusexpect(n,p,cn,1))<mtf)
       mincn <- cn
     else
       maxcn <- cn
     if (maxcn <= mincn+1)
       break
  }     # repeat
  cn
}
     

randconf <- function (n,p){
  gv <- rep(FALSE, times=n)
  m <- sample(n,p)
  gv[m] <- TRUE
  gv
}

#
# Regression FPC iteration;
# output: coef=regression coefficients, var=residual variance, g=indicator vector
rfpi <- function (indep, dep, p, gv, ca, maxit, plot) {
  cachange <- (is.na(ca))
  n <- length(gv)
  ir <- 0
  gva <- rep(FALSE, times=n)
#  print("gva") 
  change <- TRUE
  while (change) {
    ir <- ir+1
    change <- !identical(gv,gva)
    if (ir>maxit){
      change <- FALSE
      warning("Maximum iteration number exceeded.")
    }
#    oon <- on
#    on <- sum(gva)
#    orv <- rv
#    if(sum(is.na(gv))>0)
#      cat("coefs :",rc,"\nresiduals :",res,"\ngv :",gv,"\ngva :",gva,
#          "\nrv : ",rv," ca= ",ca,"\noon= ",oon," on= ",on," orv= ",orv,
#          "\nogv :",ogv,"\nrv1= ",rv1,"\nres1= ",res1,"\n")
    gva <- gv
#    print("Now reg:")
#    print(length(gv))
#    print(length(dep))
    if (p>0)
      reg <- lsfit(indep, dep, wt=gv)
    else
      reg <- lsfit(indep, dep, intercept=FALSE, wt=gv)
    rc <- coef(reg)
    res <- resid(reg)
    rv <- sum(res[as.logical(gv)]^2)/(sum(gv)-p-1)
    if (cachange){
      ca <- can(sum(gv),p)
#      cat("n= ",sum(gv)," p= ",p," ca= ",ca," change= ",change,"\n")
    }  
    gv <- res^2<=ca*rv
    coll <- FALSE
    if (p>0)
      coll <- (reg$qr$rank<ncol(as.matrix(indep))+1)
    if (coll)
      break
#    print("gv")
#    print(change)
    if (plot)
      plot(cbind(indep,dep)[,c(1,p+1)],col=1+gv,pch=1+gv)
  }
  out <- list(coef=rc, var=rv, g=gv, coll=coll, ca=ca)
  out
}

#
# fixreg main and output
#
fixreg <- function (indep=rep(1,n), dep, n=length(dep),
                    p=ncol(as.matrix(indep)),
                    ca=NA, mnc=NA, mtf=3, ir=NA, irnc=NA,
                    irprob=0.95, mncprob=0.5, maxir=20000, maxit=5*n,
                    distcut=0.85, init.group=list(), 
                    ind.storage=FALSE, countmode=100, 
                    plot=FALSE){

# Initializations, WDC

# print("fixreg")
  if (is.na(ca))
    ca <- can(n,p)
  if(is.na(ir) & is.na(irnc))
    irnc <- round(n/5)
  if(is.na(ir))
      ir <- itnumber(n,p,irnc,mtf,irprob,maxir)
  if(is.na(mnc) & ir>0)
      mnc <- minsize(n,p,ir,mtf,mncprob)
  if(is.na(mnc) & ir==0)
      mnc <- 1
  if (ir<mtf)
    mtf <- 1
  tsc <- 0          # number of too small FPCs
  ncoll <- 0        # number of iterations with collinear regressors
  gv <- rep (TRUE, times=n)
  fpc <- rfpi(indep, dep, p, gv, ca, maxit, plot)   # WDC-iteration
  if (fpc$coll){
    print("Independent variables of whole data collinear. Break forced.")
    break
  }
  imatrix <- c(sum(fpc$g)) # matrix of FPC size and intersections
  smatrix <- c(sum(fpc$g)) # matrix of FPC size and dissimilarities
  nc <- 1                  # number of FPCs
  clist <- list(fpc$coef)  # list of FPC regression coefficients
  vlist <- list(fpc$var)   # list of FPC residual variances
  nfound <- c(1)           # times found per FPC
  expectratio <- c()       # nfound/nfound expected
  if (ind.storage)
    glist <- list(fpc$g)   # list of FPC indicator vectors

# Standard iterations

  i <- 1
# print("standard it.")
  while(i<=ir){
    if(countmode*round(i/countmode)==i)
      cat("Iteration run no. ", i, " of ",ir, "\n")
    gv <- randconf(n,p+2)
# cat(length(gv)," ",sum(gv)," ")
# print("iteration")
    fpc <- rfpi(indep, dep, p, gv, ca, maxit, plot)
# print("neu?")
    neu <- TRUE
    if(!fpc$coll){
      j <- 1
      while(j<=nc) {
        if (identical(clist[[j]],fpc$coef) & identical(vlist[[j]],fpc$var)) {
          neu <- FALSE
          cnum <- j
        }          # if j found
        j <- j+1
      }            # for j
      ng <- sum(fpc$g)
      if (neu & (ng>=mnc)){
# print("neu! smatrix")
        nc <- nc +1
        nfound[nc] <- 1
        clist[[nc]] <- fpc$coef
        vlist[[nc]] <- fpc$var
# cat("ind.storage= ",ind.storage,  "nc= ", nc, "\n")
        if (ind.storage){
          glist[[nc]] <- fpc$g
          j <- 1
          while( j<=(nc-1)) {
            imatrix[sseg(nc,j)] <- sum(fpc$g & glist[[j]])
# cat("   (storage T) smatrix: ",nc,j,smatrix[sseg(nc,j)],"\n")
            j<- j+1
          }        # for j
        }          # if ind.storage
        else{
          for(j in 1:(nc-1)){
            glistj <- ((dep - cbind(1,indep) %*% clist[[j]])^2 <= ca*vlist[[j]])
            imatrix[sseg(nc,j)] <- sum(fpc$g & glistj)
# cat("   smatrix: ",nc,j,smatrix[sseg(nc,j)],"\n")
          }       # for j
        }          # else (!ind.storage)
        imatrix[sseg(nc,nc)] <- sum(fpc$g)
      }            # if neu & ng>=mnc
      if (ng<mnc)
        tsc <- tsc+1
      if (!neu)
        nfound[cnum] <- nfound[cnum]+1
    } # if !coll
    else
      ncoll <- ncoll+1
    i <- i+1
# print("schleifenende")
  }              # while i
        
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
      fpc <- rfpi(indep, dep, p, gv, ca, maxit, plot)
#      print("fpc")
      neu <- TRUE
      if(!fpc$coll){
        cnum <- nc +1
        for (j in 1:nc) {
          if (identical(clist[[j]],fpc$coef) & identical(vlist[[j]],fpc$var)) {
            neu <- FALSE
            cnum <- j
          }         # if j found
#        print(j)
        }           # for j
        ng <- sum(fpc$g)
        if (neu & (ng>=mnc)){
          nc <- nc +1
          nfound[nc] <- 1
          clist[[nc]] <- fpc$coef
          vlist[[nc]] <- fpc$var
          if (ind.storage){
            glist[[nc]] <- fpc$g
            for(j in 1:(nc-1))
              imatrix[sseg(nc,j)] <- sum(fpc$g & glist[[j]])
          }         # if ind.storage
          else{
            for(j in 1:(nc-1)){
              glistj <- ((dep - cbind(1,indep) %*% clist[[j]])^2 <= ca*vlist[[j]])
              imatrix[sseg(nc,j)] <- sum(fpc$g & glistj)
            }       # for j
          }         # else (!ind.storage)
        imatrix[sseg(nc,nc)] <- sum(fpc$g)
        }           # if neu & ng>=mnc
        if (ng<mnc)
          tsc <- tsc+1
        if(!neu)
          nfound[cnum] <- nfound[cnum]+1
      } # if !coll
      else{
        ncoll <- ncoll+1
        cnum <- 0
      } # else (coll)
      grfpc[i] <- cnum
      i <- i+1
    }             # for i
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
#  print(struc[1])
  stn <- max(struc)
  rm(comat)

# print("description of structures")

  tf <- rep(0, times=stn)   # structure: times found
  maxf <- rep(0, times=stn) # structure, rep. fpc: expectratio
  sfpc <- rep(0, times=stn) # structure: representative fpc
  for(i in 1:nc){
    expect <- clusexpect(n, p, imatrix[sseg(i,i)], ir)
    expectratio[i] <- nfound[i] / expect
    tf[struc[i]] <- tf[struc[i]] + nfound[i]
    if(expectratio[i]==maxf[struc[i]])
      if(imatrix[sseg(i,i)]<imatrix[sseg(sfpc[struc[i]],sfpc[struc[i]])])
        sfpc[struc[i]] <- i
    if(expectratio[i]>maxf[struc[i]]){
        sfpc[struc[i]] <- i
        maxf[struc[i]] <- expectratio[i]
    }
  } # for i
  stn <- sum(tf>=mtf)
#  print(stn)
  tsc <- tsc+sum(tf)-sum(tf[tf>=mtf])
  
# Output
  if (!ind.storage)
    glist <- FALSE
  out <- list(nc=nc, g=glist, coefs=clist, vars=vlist,
                nfound=nfound, er=expectratio, tsc=tsc, ncoll=ncoll,
                grto=grfpc, imatrix=imatrix,
                smatrix=smatrix, 
                stn=stn, stfound=tf,
                sfpc=sfpc, ssig=sfpc[tf>=mtf], sto=order(-tf), 
                struc=struc, n=n, p=p, ca=ca, ir=ir,
                mnc=mnc, mtf=mtf, distcut=distcut)  
  class(out) <- "rfpc"
# print ("fixreg ende")
  out   
}

summary.rfpc <- function(object, ...){
#  print("Beginn reducereg")
  clist <- list(0)
  vlist <- list(0)
  expectratio <- c()
  tf <- c()
  sn <- c()
  strn <- object$stn
#  print(strn)
  if (strn>0)
    for(i in 1:strn){
  #  cat(i, object$coefs[[object$sfpc[object$sto[i]]]], "\n")
      clist[[i]] <- object$coefs[[object$sfpc[object$sto[i]]]]
  #  cat(i, object$vars[[object$sfpc[object$sto[i]]]], "\n")
      vlist[[i]] <- object$vars[[object$sfpc[object$sto[i]]]]
  #  cat(i, object$stfound[object$sto[i]], "\n")
      tf[i] <- object$stfound[object$sto[i]]
      sn[i] <- object$imatrix[sseg(object$sfpc[object$sto[i]],object$sfpc[object$sto[i]])]
      expect <- clusexpect(object$n, object$p, sn[i], object$ir)
      expectratio[i] <- tf[i] / expect
    }
#  print(expectratio)
  sim <- simmatrix(object)
  out <- list(coefs=clist, vars=vlist, stn=strn, stfound=tf, sn=sn, 
              ser=expectratio, tsc=object$tsc, sim=sim,
              ca=object$ca, ir=object$ir, mnc=object$mnc, mtf=object$mtf)
  class(out) <- "summary.rfpc"
  out
}

fpclusters.rfpc <- function(object, indep=NA, dep=NA, ca=object$ca, ...){
  glist <- list()
  if (identical(indep,NA) & !identical(dep,NA))
    indep <- rep(1,object$n)
# cat("stn= ",object$stn,"\n")
  if(object$stn>0)
    for(i in 1:object$stn){
      if (object$g==FALSE){
        rc <- object$coefs[[ object$sfpc[object$sto[i]] ]]
        if (object$p==0)
          glist[[i]] <- (dep - rc)^2 <=
                         ca*object$vars[[ object$sfpc[object$sto[i]] ]]
        else
          glist[[i]] <- (dep - cbind(1, indep) %*% rc)^2 <=
                         ca*object$vars[[ object$sfpc[object$sto[i]] ]]
        if(sum(glist[[i]])<length(object$coefs[[1]])+1){
          cat("Warning! FPC ",i," too small, presumably because of rounding error\n")
          cat("To get the correct indicator vector, run fixreg again with 'ind.storage=T'\n") 
        } # if sum glist too small
      }   # if g==F
      else
        glist[[i]] <- object$g[[object$sfpc[object$sto[i]]]]
  # print(i)
    } # for i
  else
    warning("No FPCs were found often enough.")
  glist
}

# Visualization of the representative FPC no. no; x-axis: regression 
# linera combination => theoretical line= identity.
plot.rfpc <- function(x, indep=rep(1,n), dep, no, bw=TRUE,
                      main=c("Representative FPC No. ",no),
                      xlab="Linear combination of independents",
                      ylab=deparse(substitute(indep)),
                      xlim=NULL, ylim=range(dep), 
                      pch=NULL, col=NULL,...){
  n <- x$n
  p <- x$p
  sumobj <- summary(x)
  rc <- sumobj$coefs[[no]]
  lc <- cbind(1, indep) %*% rc
  ind <- ((dep-lc)^2 <= sumobj$ca*sumobj$vars[[no]])
  if (is.null(pch))
    pch <- ifelse(bw, 2, 1)
  if (is.null(col))
    col <- ifelse(bw, 1, 2)
  if (is.null(xlim))
    xlim <- range(lc)
  if (is.null(ylim))
    ylim <- range(dep)
  plot(lc[ind], dep[ind], main=main, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim, pch=pch, col=col, ...)
  points(lc[!ind], dep[!ind], pch=1, col=1, ...)
  abline(c(0,1))
  abline(c(-sqrt(sumobj$ca*sumobj$vars[[no]]),1), lty="dotted")
  abline(c(sqrt(sumobj$ca*sumobj$vars[[no]]),1), lty="dotted")
  invisible()
}

print.rfpc <- function(x, ...){
  cat("Linear Regression Fixed Point Cluster object\n")  
  cat(x$stn," representative stable fixed point clusters\n") 
  cat(" of totally ",x$nc," found fixed point clusters.\n")
  invisible(x)
}

# Summary output of rfpc objects; choose fpcobj <- summary(fpcobj) to
# replace fpcobj by its reduction    
print.summary.rfpc <- function(x, maxnc=30, ...){
#  print("Beginn Summary")
  minnc <- min(x$stn,maxnc)
  cat("  *  Fixed Point Clusters  *\n\n")
  cat("Often a clear cluster in the data leads to several similar FPCs.\n")
  cat("The summary shows the representative FPCs of groups of similar FPCs,\n")
  cat("which were found at least ",x$mtf," times.\n\n")
  cat("Constant ca= ",x$ca,"\n")
  cat("Number of representative FPCs: ", x$stn, "\n\n")
  cat("FPCs with less than ",x$mnc," points were skipped.\n")
  cat(x$tsc," iterations led to skipped FPCs.\n\n")
  if (x$stn>maxnc)
    cat("Warning! Only ",maxnc," clusters are displayed. \nSpecify maxnc in call of summary.rfpc if you want enlarged display or \nrun fixreg with larger ca, mnc or mtf to get fewer clusters. \n\n")
  if (x$stn==0)
    cat("No FPCs were found often enough.\n")
  else{
    for(i in 1:minnc){
      cat(" FPC ",i, "\n")
      cat("  Times found (group members): ",x$stfound[i], "\n")
      if (x$ir>0)
        cat("  Ratio to estimated expectation: ",x$ser[i], "\n")
      cat("  Regression parameters:\n")
      print(x$coefs[[i]])
      cat("  Error variance: ", x$vars[[i]], "\n")
      cat("  Number of points: ",x$sn[[i]], "\n\n")
    }
    cat("Number of points in intersection of  representative FPCs\n")
    sm <- rep(0,times=minnc^2)
    dim(sm) <- c(minnc,minnc)
    for(i in 1:minnc)
      for(j in 1:minnc)
        sm[i,j] <- x$sim[sseg(i,j)]
    print(sm)
  }
  invisible(x)
}

