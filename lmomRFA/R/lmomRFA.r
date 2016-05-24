#***********************************************************************
#*                                                                     *
#*  R code written for the R package "lmomRFA"                         *
#*                                                                     *
#*  J. R. M. HOSKING <jrmhosking@gmail.com>                            *
#*                                                                     *
#*  Version 3.0-1  February 2015                                       *
#*                                                                     *
#***********************************************************************

cluagg<-function(x, method="ward") {
  dx<-dist(x)
  if (pmatch(method,"ward",nomatch=0)==1) dx<-dx^2
  hc<-hclust(dx,method=method)

  hcm<-(-hc$merge)
  for (j in 1:nrow(hcm)) {
    if (hcm[j,1]<0) hcm[j,1]<-hcm[-hcm[j,1],1]
    if (hcm[j,2]<0) hcm[j,2]<-hcm[-hcm[j,2],1]
    hcm[j,]<-sort(hcm[j,])
  }
  return(list(merge=hcm,wgss=cumsum(hc$height/2)))
}

cluinf<-function(merge, nclust) {
  if (is.list(merge) && names(merge)==c("merge","wgss")) merge<-merge$merge
  n<-nrow(merge)+1
  vec<-1:n
  for (i in 1:(n-nclust)) vec[vec==merge[i,2] ] <- vec[merge[i,1] ]
  assign<-match(vec,sort(unique(vec)))
  list<-lapply(1:nclust,function(i) which(assign==i))
  padnum<-formatC(1:nclust,width=nchar(nclust),format="d",flag="0")
  names(list)<-paste("cluster.",padnum,sep="")
  num<-sapply(list,length)
  return(list(assign=assign,list=list,num=num))
}

clukm<-function(x, assign, maxit=10, algorithm="Hartigan-Wong") {
  x<-as.matrix(x)
  if (nrow(x)!=length(assign))
    stop("number of rows of 'x' and length of 'assign' must be equal")
  centers<-apply(x,2,function(y) tapply(y,assign,mean))
  stats::kmeans(x,centers,iter.max=maxit,algorithm=algorithm)
}

reglmr<-function(xmom, weight) {
## Regional weighted average of L-moments
  xmom<-as.matrix(xmom)
  if (missing(weight)) weight<-rep(1,nrow(xmom))
  if (length(weight)!=nrow(xmom))
    stop("number of rows of 'xmom' and length of 'weight' must be equal")
  if (ncol(xmom)>1) xmom[,2]<-xmom[,2]/xmom[,1]
  xmom[,1]<-1
  apply(xmom,2,weighted.mean,w=weight,na.rm=TRUE)
}

as.regdata<-function(x, warn.names=TRUE) {
##  Convert an R object to class "regdata"
  if (inherits(x,"regdata")) return(x)
#
  x<-as.data.frame(x,row.names=NULL)
# Sanity checks for the data in 'x'
  if (ncol(x)<=4) stop("must have at least 4 columns")
  sitenames<-x[,1]<-as.character(x[,1])
  if (any(duplicated(sitenames[!is.na(sitenames)])))
    warning("site names are not all different")
  if (any(x[,2]!=round(x[[2]]))) stop("record lengths are not all integers")
  if (any(x[,2]<=0)) stop("record lengths are not all greater than zero")
  if (any(x[,3]<=0)) stop("site means are not all greater than zero")
  if (any(x[,4]<0 | x[,4]>1))
    stop("L-CV values are not all between 0 and 1")
  if (any(x[,-(1:4)]<(-1) | x[,-(1:4)]>1))
    stop("L-moment ratios are not all between -1 and +1")
# Check names of columns of 'x'
  nx<-names(x)
  if (identical(nx,paste("V",seq(along=x),sep="")))
    names(x)<-c("name","n","mean","t",paste("t_",seq(3,length(x)-2),sep=""))
  else if (warn.names) {
    mygrep<-function(x,y) (regexpr(x,y,ignore.case=TRUE)>0)
    if (!mygrep("^(id|name|site|site[[:punct:]]?id|site[[:punct:]]?name)s?$",nx[1]))
      warning('Column 1 (site name) has nonstandard name "',nx[1],'"')
#   if (!mygrep("",nx[2]))
#     warning('Column 2 (record length) has nonstandard name "',nx[2],'"')
    if (!mygrep("^(mean|l[[:punct:]]?1)$",nx[3]))
      warning('Column 3 (mean) has nonstandard name "',nx[3],'"')
    if (!mygrep("^(t|t[[:punct:]]2|lcv|l[[:punct:]]?cv)$",nx[4]))
      warning('Column 4 (L-CV) has nonstandard name "',nx[4],'"')
    if (length(x)>=5 && !mygrep("(^l[[:punct:]]s)|(^t[[:punct:]]?3$)$",nx[5]))
      warning('Column 5 (L-skewness) has nonstandard name "',nx[5],'"')
    if (length(x)>=6 && !mygrep("(^l[[:punct:]]k)|(^t[[:punct:]]?4$)$",nx[6]))
      warning('Column 6 (L-kurtosis) has nonstandard name "',nx[6],'"')
    if (length(x)>=7)
      for (j in 7:length(x)) {
        if(!mygrep(paste("^t[[:punct:]]?",j-2,"$",sep=""),nx[j]))
        warning('Column ',j,' (t-',j-2,') has nonstandard name "',nx[j],'"')
      }
  }
  return(structure(x,class=c("regdata",class(x))))
}

regsamlmu<-function(x, nmom=5, sort.data=TRUE, lcv=TRUE) {
  if (is.list(x)) {
    xmom<-t(sapply(x,samlmu,nmom=nmom,sort.data=sort.data))
    n<-sapply(x,function(y) length(y[!is.na(y)]))
    name<-names(x)
    if (is.null(name)) name<-seq_along(x)
  } else {
    x<-as.matrix(x)
    xmom<-t(apply(x,2,samlmu,nmom=nmom,sort.data=sort.data))
    n<-apply(x,2,function(y) length(y[!is.na(y)]))
    name<-colnames(x)
    if (is.null(name)) name<-seq_len(ncol(x))
  }
  if (nmom==0) xmom<-matrix(nrow=length(name),ncol=0)
  if (nmom==1) xmom<-matrix(xmom,ncol=1,dimnames=list(NULL,"l_1"))
  if (lcv && nmom>=2) {
    xmom[,2]<-xmom[,2]/xmom[,1]
    colnames(xmom)[2]<-"t"
  }
  out<-cbind(data.frame(name=name,stringsAsFactors=FALSE),n,xmom)
  row.names(out)<-NULL
  if (sort.data && lcv && nmom>=2) class(out)<-c("regdata",class(out))
  return(out)
}

regavlmom<-function(regdata, weight) {
## Regional weighted average of L-moments
  regdata<-as.regdata(regdata)
  if (missing(weight)) weight<-regdata[[2]]
  if (length(weight)==1) weight<-rep(weight,length=nrow(regdata))
  if (length(weight)!=nrow(regdata))
    stop("number of rows of 'regdata' and length of 'weight' must be equal")
  rmom<-apply(regdata[-(1:3)],2,weighted.mean,w=weight,na.rm=TRUE)
  rmom<-c(1,rmom)
  names(rmom)[1:2]<-c("l_1","l_2")
  return(rmom)
}

regtst<-function(regdata, nsim=1000){
##  Discordancy, heterogeneity and goodness-of-fit statistics
##  for regional frequency analysis

  regdata<-as.regdata(regdata)

  if (ncol(regdata)<6)
    stop("'regdata' must be a data frame with at least 6 columns")

  nsites<-nrow(regdata)
  len<-regdata[,2]
  maxrec<-max(len)
  xmom <- if (ncol(regdata)==6) t(cbind(regdata[,3:6],2)) else t(regdata[,3:7])

  fort<-.Fortran("regtst",PACKAGE="lmomRFA",
    nsites=as.integer(nsites),
    len=as.integer(len),
    xmom=as.double(xmom),
    nsim=as.integer(nsim),
    rmom=double(5),
    d=double(nsites),
    vobs=double(3),
    vbar=double(3),
    vsd=double(3),
    h=double(3),
    z=double(5),
    para=double(30),
    rpara=double(4),
    t4fit=double(5),
    work=double(nsites*3),
    x=double(maxrec),
    maxrec=as.integer(maxrec))

  if (all(fort$d==0)) {
    is.na(fort$d[])<-TRUE
    warning("unable to invert sum-of-squares matrix - D statistics not calculated")
  }

  if (nsim<=1) fort[c("rpara","vobs","vbar","vsd","h","z","t4fit")]<-NULL

  dc1<-c(3,3,3,3,1.3330,1.6481,1.9166,2.1401,2.3287,2.4906,
    2.6321,2.7573,2.8694,2.9709,3,3,3,3)
  dc2<-c(4,4,4,4,1.3333,1.6648,1.9821,2.2728,2.5337,2.7666,
    2.9748,3.1620,3.3310,3.4844,3.6246,3.7532,3.8718,3.9816)
  Dcrit <- if (nsites>length(dc1)) c(3,4) else c(dc1[nsites],dc2[nsites])

  para=list(
    glo=fort$para[1:3],
    gev=fort$para[6:8],
    gno=fort$para[11:13],
    pe3=fort$para[16:18],
    gpa=fort$para[21:23],
    wak=fort$para[26:30])
  names(para$glo)<-lmom:::lmom.dist$glo$parnames
  names(para$gev)<-lmom:::lmom.dist$gev$parnames
  names(para$gno)<-lmom:::lmom.dist$gno$parnames
  names(para$pe3)<-lmom:::lmom.dist$pe3$parnames
  names(para$gpa)<-lmom:::lmom.dist$gpa$parnames
  names(para$wak)<-lmom:::lmom.dist$wak$parnames

  if (ncol(regdata)==6) {fort$rmom[5]<-NA; para$wak[]<-NA}

  out<-list(
    data=regdata,
    nsim=nsim,
    D=fort$d,
    Dcrit=Dcrit,
    rmom=fort$rmom,
    rpara=fort$rpara,
    vobs=fort$vobs,
    vbar=fort$vbar,
    vsd=fort$vsd,
    H=fort$h,
    para=para,
    t4fit=fort$t4fit,
    Z=fort$z)

  names(out$D)<-regdata[[1]]
  names(out$rmom)<-c("mean","t","t_3","t_4","t_5")
  if (nsim>1) {
    names(out$rpara)<-lmom:::lmom.dist$kap$parnames
    names(out$t4fit)<-names(out$Z)<-c("glo","gev","gno","pe3","gpa")
  }

  class(out)<-"regtst"
  return(out)
}

regtst.s<-function(regdata, nsim=1000) {
##
##  (mostly) native S version of regtst()
##
##  "Mostly", because calculations in inner loop use fast versions
##  of quakap() and samlmu() that call Fortran routines.
##
  regdata<-as.regdata(regdata)

  if (ncol(regdata)<6)
    stop("'regdata' must be a data frame with at least 6 columns")

  nsites<-nrow(regdata)
  len<-regdata[,2]
  xmom <- if (ncol(regdata)==6) cbind(as.matrix(regdata[,4:6]),NA) else as.matrix(regdata[,4:7])

  # Compute discordancy measure
  if (nsites<=3) D<-rep(1,nsites)
  else {
    u<-xmom[,1:3]
    D<-try(nsites/(nsites-1)/3*mahalanobis(u,colMeans(u),cov(u)),silent=TRUE)
    if (inherits(D,"try-error")) {
      D<-rep(NA_real_,nsites)
      warning("unable to invert sum-of-squares matrix - D statistics not calculated")
    }
  }
  names(D)<-regdata[[1]]

  dc1<-c(3, 3, 3, 3, 1.333, 1.6481, 1.9166, 2.1401, 2.3287,
    2.4906, 2.6321, 2.7573, 2.8694, 2.9709, 3, 3, 3, 3)
  dc2<-c(4, 4, 4, 4, 1.3333, 1.6648, 1.9821, 2.2728, 2.5337,
    2.7666, 2.9748, 3.162, 3.331, 3.4844, 3.6246, 3.7532,
    3.8718, 3.9816)
  Dcrit <- if (nsites>length(dc1)) c(3,4) else c(dc1[nsites],dc2[nsites])

  rmom<-c(1,apply(xmom,2,weighted.mean,w=len))

  para<-list(glo=pelglo(rmom), gev=pelgev(rmom), gno=pelgno(rmom),
             pe3=pelpe3(rmom), gpa=pelgpa(rmom),
             wak = if (ncol(regdata)==6)
               setNames(rep(NA_real_,5),lmom:::lmom.dist$wak$parnames) else pelwak(rmom))
# or... #    wak = if (ncol(regdata)==6) {z<-pelwak(c(0,1,0,0,0));is.na(z)<-TRUE;z} else pelwak(rmom))

  if (nsim<=1) {

    rpara<-vobs<-vbar<-vsd<-H<-Z<-t4fit<-NULL

  } else {

    rpara <- if (rmom[4]<=(1+5*rmom[3]^2)/6) pelkap(rmom) else c(pelglo(rmom),-1)

    # Fast version of quakap(), for use in inner loop
    my.quakap<-function(f,para) .Fortran("qkap",PACKAGE="lmomRFA",
        as.double(f),length(f),as.double(para))[[1]]

    re<-replicate(nsim, {
      lmrsim<-sapply(seq_len(nsites), function(j) {  # For each site:
        nrec<-len[j]
        simdat<-my.quakap(runif(nrec),rpara)         # - Generate simulated data
        .samlmu(simdat)                              # - Compute L-moment ratios
      })
      lmrsim[2,]<-lmrsim[2,]/lmrsim[1,]              # Compute L-CV of simulated data
      regave<-apply(lmrsim,1,weighted.mean,w=len)
      sw<-sweep(lmrsim,1,regave)
      v1<-sqrt(weighted.mean(sw[2,]^2,w=len))
      v2<-weighted.mean(sqrt(sw[2,]^2+sw[3,]^2),w=len)
      v3<-weighted.mean(sqrt(sw[3,]^2+sw[4,]^2),w=len)
      c(v1,v2,v3,regave[4])    # V statistics and regional average t_4
    })
    re.bar<-apply(re,1,mean)
    re.sd<-apply(re,1,sd)
    vbar<-re.bar[1:3]
    vsd<-re.sd[1:3]
    t4bias<-re.bar[4]-rmom[4]
    t4sd<-re.sd[4]

    sw<-sweep(t(xmom[,1:3]),1,rmom[2:4])
    v1<-sqrt(weighted.mean(sw[1,]^2,w=len))
    v2<-weighted.mean(sqrt(sw[1,]^2+sw[2,]^2),w=len)
    v3<-weighted.mean(sqrt(sw[2,]^2+sw[3,]^2),w=len)
    vobs<-c(v1,v2,v3)
    H<-(vobs-vbar)/vsd

    t4fit<-sapply(names(para)[1:5],function(dist) lmom:::lmrxxx(dist,para[[dist]],4)$lmom[4])
    Z<-(t4fit-rmom[4]+t4bias)/t4sd

    names(rpara)<-lmom:::lmom.dist$kap$parnames
    names(t4fit)<-names(Z)<-names(para)[1:5]
  }

  out<-list(data=regdata, nsim=nsim, D=D, Dcrit=Dcrit,
    rmom=rmom, rpara=rpara, vobs=vobs, vbar=vbar, vsd=vsd, H=H,
    para=para, t4fit=t4fit, Z=Z)
  names(out$rmom)<-c("mean", "t", "t_3", "t_4", "t_5")
  class(out)<-"regtst"
  return(out)
}

print.regtst<-function(x,...) {
## Print method for an object of class "regtst"
  cat("Discordancy measures (critical value ",formatC(x$Dcrit[1],2,format="f"),")\n",sep="")
  cat(formatC(x$D,2,format="f"),"\n\n")
  if (x$nsim<=1) {
    cat("Heterogeneity measures not calculated\n")
    cat("Goodness-of-fit measures not calculated\n")
  } else {
    cat("Heterogeneity measures (based on",x$nsim,"simulations)\n")
    cat(formatC(x$H,2,format="f"),"\n\n")
    cat("Goodness-of-fit measures (based on",x$nsim,"simulations)\n")
    print(round(x$Z,2))
  }
  return(invisible(x))
}

summary.regtst<-function(object,
  prob=c(0.01,0.02,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.98,0.99,0.999),
  conf=0.90, decimals=c(4,4,2,3), ...){
## Summary method for an object of class "regtst"

  prob<-prob[prob>=0 & prob<=1]

  if (is.null(object$para)) {
    out<-c(object,list(prob=prob,quant=NULL,decimals=decimals))
  } else {
    quant<-matrix(NA_real_,nrow=6,ncol=length(prob))
    if (ncol(quant)>0) {
      quant[1,]<-quaglo(prob,object$para$glo)
      quant[2,]<-quagev(prob,object$para$gev)
      quant[3,]<-quagno(prob,object$para$gno)
      quant[4,]<-quape3(prob,object$para$pe3)
      quant[5,]<-quagpa(prob,object$para$gpa)
      if (!is.na(object$para$wak[1])) quant[6,]<-quawak(prob,object$para$wak)
      colnames(quant)<-format(prob,scientific=FALSE)
    }
    rownames(quant)<-names(object$para)
    out<-c(object,list(conf=conf,prob=prob,quant=quant,decimals=decimals))
  }

  class(out)<-"summary.regtst"
  return(out)

}

print.summary.regtst<-function(x, decimals, ...) {
## Print an object of class "summary.regtst"
  if (missing(decimals)) decimals<-x$decimals
  ldec<-length(decimals)
  if (ldec<4) decimals<-c(decimals,c(4,4,2,3)[(ldec+1):4])
  dlmom<-decimals[1]
  dpara<-decimals[2]
  dtest<-decimals[3]
  dquant<-decimals[4]

  dat<-x$data
  dat[,-(1:2)]<-format(dat[,-(1:2)],digits=1,nsmall=dlmom,scientific=FALSE)
  nsites<-length(x$D)
  stars<-rep("  ",nsites)
  substring(stars,1,1)[x$D>=x$Dcrit[1] ]<-"*"
  substring(stars,2,2)[x$D>=x$Dcrit[2] ]<-"*"
  print(cbind(dat,"D(i)"=format(x$D,digits=1,nsmall=dtest),"  "=stars),
    right=FALSE)

  cat("\n")
  rmommat<-matrix(x$rmom[-1],nrow=1,
    dimnames=list("Weighted means  ",names(x$rmom)[-1]))
  print(noquote(formatC(rmommat,digits=dlmom,format="f")))

  if (any(stars!="  "))
    cat("\nFlagged test values:",
      formatC(sort(x$D[stars!="  "],decreasing=TRUE),digits=2,format="f"))

  if (x$nsim>1) {
    cat("\n")
    parmat<-matrix(x$rpara,nrow=1,
      dimnames=list("Parameters of regional kappa distribution  ",names(x$rpara)))
    print(noquote(formatC(parmat,digits=dpara,format="f")))
  }

  dnames<-c(
    "Gen. logistic      ",
    "Gen. extreme value ",
    "Gen. normal        ",
    "Pearson type III   ",
    "Gen. Pareto        ",
    "Wakeby             ")

  conf<-NA

  if (x$nsim>1) {

    if (nsites>1) {
      vobs<-formatC(x$vobs,digits=dlmom,width=dlmom+3,format="f")
      vbar<-formatC(x$vbar,digits=dlmom,width=dlmom+3,format="f")
      vsd <-formatC(x$vsd ,digits=dlmom,width=dlmom+3,format="f")
      H   <-formatC(x$H   ,digits=dtest,width=dtest+3,format="f")
      stars<-paste(format("",width=max(0,dlmom-dtest)),
        ifelse(x$H>2,"**",ifelse(x$H>1,"* ","  ")),sep="")
      cat("\n\n*****  HETEROGENEITY MEASURES  *****\n")
      cat("Number of simulations =",x$nsim,"\n\n")
      cat("Observed     s.d. of group L-CV             =",vobs[1],"\n")
      cat("Sim. mean of s.d. of group L-CV             =",vbar[1],"\n")
      cat("Sim. s.d. of s.d. of group L-CV             =",vsd[1] ,"\n")
      cat("Heterogeneity measure H[1]                  =",H[1],stars[1],"\n\n")
      cat("Observed     s.d. of L-CV / L-skew distance =",vobs[2],"\n")
      cat("Sim. mean of s.d. of L-CV / L-skew distance =",vbar[2],"\n")
      cat("Sim. s.d. of s.d. of L-CV / L-skew distance =",vsd[2] ,"\n")
      cat("Heterogeneity measure H[2]                  =",H[2],stars[2],"\n\n")
      cat("Observed     s.d. of L-skew/L-kurt distance =",vobs[3],"\n")
      cat("Sim. mean of s.d. of L-skew/L-kurt distance =",vbar[3],"\n")
      cat("Sim. s.d. of s.d. of L-skew/L-kurt distance =",vsd[3] ,"\n")
      cat("Heterogeneity measure H[3]                  =",H[3],stars[3],"\n\n")
    }

    stars<-rep(" ",5)
    if (is.numeric(x$conf) && length(x$conf)==1 && x$conf>0 && x$conf<1) {
      conf<-format(x$conf,nsmall=2)
      Zcrit<-qnorm(0.5+x$conf/2)
      stars[abs(x$Z)<=Zcrit]<-"*"
    }

    t4fit<-formatC(x$t4fit,digits=dlmom,width=dlmom+3,format="f")
    Z    <-formatC(x$Z    ,digits=dtest,width=dtest+4,format="f")
    cat("\n*****  GOODNESS-OF-FIT MEASURES  *****\n")
    cat("Number of simulations =",x$nsim,"\n\n")
    for (j in 1:5)
      cat(dnames[j],"  L-kurtosis =",t4fit[j],"   Z value =",Z[j],stars[j],"\n")
    cat("\n")

  }

  if (is.na(conf)) {
    cat("\nPARAMETER ESTIMATES\n\n")
    ok<-1:6
  } else {
    cat("\nPARAMETER ESTIMATES FOR DISTRIBUTIONS ACCEPTED AT THE",
      conf,"LEVEL\n\n")
    ok<-c(which(stars=="*"),6)
  }

  for (j in ok) { cat(dnames[j],
    formatC(x$para[[j]],digits=dpara,width=dpara+4,format="f"),"\n")
  }

  if (ncol(x$quant)>0) {
    cat("\nQUANTILE ESTIMATES\n")
    quant<-x$quant
    rownames(quant)<-dnames
    colnames(quant)<-paste(" ",colnames(quant))
    quant<-formatC(quant,digits=dquant,width=dquant+4,format="f")
    print(noquote(quant[ok,,drop=FALSE]))
  }

  return(invisible(x))
}

regfit<-function(regdata, dist) {
## Fit a regional frequency distribution
  regdata<-as.regdata(regdata)

  if (!is.character(dist) || length(dist)!=1)
    stop("'dist' must be a character string")
  pelname<-paste("pel",dist,sep="")
  quaname<-paste("qua",dist,sep="")
  pf<-parent.frame()
  if (!exists(pelname,mode="function",envir=pf)) stop('function "',pelname,'" not found')
  if (!exists(quaname,mode="function",envir=pf)) stop('function "',quaname,'" not found')
  pelfun<-get(pelname,mode="function",envir=pf)
  quafun<-get(quaname,mode="function",envir=pf)

  rmom<-regavlmom(regdata)
  para<-pelfun(rmom)
  out<-structure(
    list(
      dist=dist,
      para=para,
      qfunc=function(f) quafun(f,para),
      rmom=rmom,
      index=structure(regdata[[3]],names=regdata[[1]])),
    class="rfd")
  return(out)
}

print.rfd<-function(x,...) {
  cat("Regional frequency distribution:",x$dist,"\nParameters:\n")
  print(x$para)
  return(invisible(x))
}

regqfunc<-function(rfd) {
## Regional growth curve (quantile function of a regional frequency distribution)
  if (!inherits(rfd,"rfd")) stop("'rfd' must be an object of class \"rfd\"")
  rfd$qfunc
}

regquant<-function(f,rfd) {
## Quantiles of the regional growth curve
  if (!inherits(rfd,"rfd")) stop("'rfd' must be an object of class \"rfd\"")
  structure(rfd$qfunc(f),names=f)
}

siteqfunc<-function(rfd, sitenames, index) {
## Quantile functions for individual sites
  if (!inherits(rfd,"rfd")) stop("'rfd' must be an object of class \"rfd\"")
  if (missing(index)) {
    index <- if (missing(sitenames)) rfd$index else rfd$index[sitenames]
  } else {
    if (!all(index>=0)) stop("values of 'index' must not be negative")
    if (!missing(sitenames)) {
      if (!is.character(sitenames))
        stop("'sitenames' must be of type \"character\" when 'index' is present")
      if (length(sitenames)!=length(index))
        stop("'sitenames' and 'index' must have the same length")
      names(index)<-sitenames
    }
  }
  out<-lapply(index,
    function(ind) eval(substitute(function(f) a*rfd$qfunc(f)),list(a=ind)))
  if (length(out)==1) out<-out[[1]]
  out
}

sitequant<-function(f, rfd, sitenames, index, drop=TRUE) {
## Quantiles for individual sites
  if (!inherits(rfd,"rfd")) stop("'rfd' must be an object of class \"rfd\"")
  if (missing(index)) {
    index <- if (missing(sitenames)) rfd$index else rfd$index[sitenames]
  } else {
    if (!all(index>=0)) stop("values of 'index' must not be negative")
    if (!missing(sitenames)) {
      if (!is.character(sitenames))
        stop("'sitenames' must be of type \"character\" when 'index' is present")
      if (length(sitenames)!=length(index))
        stop("'sitenames' and 'index' must have the same length")
      names(index)<-sitenames
    }
  }
  outer(index,regquant(f,rfd))[,,drop=drop]
}

regsimh <- function(qfunc, para, cor=0, nrec, nrep=500, nsim=500) {
## Simulate H and Z
   nsites<-length(nrec)
   nmax<-max(nrec)
   if (any(nrec<=0)) stop("record lengths must all be positive")

   if (!is.list(qfunc)) qfunc<-list(qfunc)
   if (!all(sapply(qfunc,is.function)))
     stop("'qfunc' must be a function or a list of functions")
   if (!is.element(length(qfunc),c(1,nsites)))
     stop("list 'qfunc' must have either 1 or 'length(nrec)' components")

   # Function q2qua() converts a "base R"-type quantile function (one argument for
   # each parameter of the distribution) into a "package lmom"-type quantile function
   # (two arguments, the second one containing all the parameters of the distribution)

   q2qua <- function(f) function(u,p) do.call(f,c(list(u),as.list(p)))

   # We don't use the reverse transformation, but here it is anyway
   #   qua2q <- function(f) function(u,...) mapply(f,u,mapply(c,...,SIMPLIFY=FALSE))
   # or (much faster if each arg in '...' is a single number)
   #   qua2q <- function(f) function(u,...) {
   #     if (length(c(...))==length(list(...))) f(u,c(...))
   #     else mapply(f,u,mapply(c,...,SIMPLIFY=FALSE))
   #   }

   # Convert all elements of the 'qfunc' list to "package-lmom" type quantile functions

   qfunc<-lapply(qfunc,function(f) if (length(formals(f))==2) f else q2qua(f))

   # qflocal() is a function with 3 arguments that generates a function call
   # its first argument with

   qflocal <- function(f,u,p) f(u,p)

   # Convert the 'para' argument into a list of (either 1 or) 'nsites' vectors;
   # each vector contains the distribution parameters for one site (or all sites).

   if (missing(para)) para<-list(NULL)
   else if (is.vector(para) && is.numeric(para)) para<-list(para)
   else {
     if (is.data.frame(para)) para<-as.matrix(para)
     if (is.matrix(para)) {
       if (!is.element(nrow(para),c(1,nsites)))
         stop("matrix or data frame 'para' must have either 1 or 'length(nrec)' rows")
       para<-lapply(split(para,row(para)),function(x) x[!is.na(x)])
     } else if (is.list(para)) {
       if (!is.element(length(para),c(1,nsites)))
         stop("list 'para' must have either 1 or 'length(nrec)' components")
       para<-lapply(para,unlist)
     } else stop("'para' must be a vector, matrix, data frame or list")
   }

   # Generate the correlation matrix if necessary, and find its Cholesky decomposition

   nocorr<-identical(cor,0)
   if (!nocorr) {
     if (is.matrix(cor)) cor<-cov2cor(cor)
       else {
         avcor<-as.vector(cor)
         if (length(avcor)!=1) stop("'cor' must be either a matrix or a scalar")
        cor<-diag(1-avcor,nsites)+avcor
      }
      cholcor<-try(chol(cor),silent=TRUE)
      if (class(cholcor)=="try-error") stop("Correlation matrix is not positive definite")
    }

   # Simulation loop

   re<-replicate(nrep, {

     # Generate uniform random variates at each site

     if (nocorr) {

       # If there is no inter-site correlation, generate independent
       # uniform samples at each site

       ulist<-lapply(nrec,runif)

     }  else {

       # If there is correlation, generate correlated normal samples
       # each of length 'nmax' ...

       zmat<-matrix(rnorm(nsites*nmax),nsites,nmax)
       zmat<-crossprod(cholcor,zmat)

       # ... transform to uniform ...

       zmat<-pnorm(zmat)

       # ... and retain only the required number of data values at each site

       ulist<-lapply(1:nsites, function(j) zmat[j,1:nrec[j] ])
     }

     # Transform the uniform variates to the required parent distribution.
     # At site i, feed the uniform sample ulist[[i]] through quantile
     # function qfunc[[i]] with parameters para[[i]].

     datlist<-mapply(qflocal,qfunc,ulist,para,SIMPLIFY=FALSE)

     # Compute the test statistics for the simulated region's data

     rt<-regtst(regsamlmu(datlist),nsim=nsim)

     # Extract the H and Z measures from the output of regtst()

     c(H=rt$H,Z=rt$Z)
   })

   # Compute averages, across simulations, of H and Z measures

   means<-rowMeans(re)

   # Return the results in an object of class "regsimh"

   return(structure(list(nrep=nrep,nsim=nsim,results=re,means=means),class="regsimh"))

}

print.regsimh<-function(x,...) {
  # Print an object of class "regsimh"
  # Print just the average H and Z measures, to 2 decimal places
  cat("Average heterogeneity measures (based on", x$nsim, "simulations within each of", x$nrep, "simulated regions)\n")
  cat(formatC(x$means[1:3], 2, format="f"), "\n\n")
  cat("Average goodness-of-fit measures (based on", x$nsim, "simulations within each of", x$nrep, "simulated regions)\n")
  names(x$means)[4:8]<-c("glo","gev","gno","pe3","gpa")
  print(round(x$means[4:8], 2))
}

regsimq<-function(qfunc, para, cor=0, index=NULL, nrec, nrep=10000,
  fit="gev", f=c(0.01,0.1,0.5,0.9,0.99,0.999), boundprob=c(0.05,0.95),
  save=TRUE) {
## Simulations for error bounds of regional growth curve

  nsites<-length(nrec)
  nmax<-max(nrec)
  if (any(nrec<=0)) stop("record lengths must all be positive")

  if (!is.list(qfunc)) qfunc<-list(qfunc)
  if (!all(sapply(qfunc,is.function)))
    stop("'qfunc' must be a function or a list of functions")
  if (!is.element(length(qfunc),c(1,nsites)))
    stop("list 'qfunc' must have either 1 or 'length(nrec)' components")

  # Function q2qua() converts a "base R"-type quantile function (one argument for
  # each parameter of the distribution) into a "package lmom"-type quantile function
  # (two arguments, the second one containing all the parameters of the distribution)

  q2qua <- function(f) function(u,p) do.call(f,c(list(u),as.list(p)))

  # We don't use the reverse transformation, but here it is anyway
  # qua2q <- function(f) function(u,...) mapply(f,u,mapply(c,...,SIMPLIFY=FALSE))

  # Convert all elements of the 'qfunc' list to "package-lmom" type quantile functions

  qfunc<-lapply(qfunc, function(func)
    if (length(formals(func))==2) func else q2qua(func))

  # Make 'qfunc' a list of 'nsites' functions,
  # the quantile functions for each site.

  if (length(qfunc)==1) qfunc<-rep(qfunc,nsites)

  # Convert the 'para' argument into a list of 'nsites' vectors,
  # each containing the distribution parameters for one site.

  if (missing(para)) para<-list(NULL)
  else if (is.vector(para) && is.numeric(para)) para<-list(para)
  else {
    if (is.data.frame(para)) para<-as.matrix(para)
    if (is.matrix(para)) {
      if (!is.element(nrow(para),c(1,nsites)))
        stop("matrix or data frame 'para' must have either 1 or 'length(nrec)' rows")
      para<-lapply(split(para,row(para)),function(x) x[!is.na(x)])
    } else if (is.list(para)) {
      if (!is.element(length(para),c(1,nsites)))
        stop("list 'para' must have either 1 or 'length(nrec)' components")
      para<-lapply(para,unlist)
    } else stop("'para' must be a vector, matrix, data frame or list")
  }
  if (length(para)==1) para<-rep(para,nsites)

  # Generate the correlation matrix if necessary, and find its Cholesky decomposition

  nocorr<-identical(cor,0)
  if (!nocorr) {
    if (is.matrix(cor)) {
      cor<-cov2cor(cor)
      if (!all(dim(cor)==nsites)) stop("matrix 'cor' must be square and of order 'length(nrec)'")
    } else {
      avcor<-as.vector(cor)
      if (length(avcor)!=1) stop("'cor' must be either a matrix or a scalar")
      cor<-diag(1-avcor,nsites)+avcor
    }
    cholcor<-try(chol(cor),silent=TRUE)
    if (class(cholcor)=="try-error") stop("Correlation matrix is not positive definite")
  }

  # Compute the index flood values if necessary

  if (is.null(index)) {
    integrator<-function(f,p) integrate(f,lower=0,upper=1,p)$value
    index<-tryCatch(mapply(integrator,qfunc,para), error=function(...) NA)
    if (any(is.na(index)))
      stop("integration failed: unable to compute index flood value at site(s) ",
        which(is.na(index)))
  } else {
   index<-as.vector(index)
   if (!is.element(length(index),c(1,nsites)))
     stop("vector 'index' must have length 1 or 'length(nrec)'")
   if (length(index)==1) index<-rep(index,length=nsites)
  }

  # Are all index flood values equal to 1?

  index1<-isTRUE(all.equal(index,rep(1,length(index))))

  # Check that the 'fit' routines exist

  if (!is.character(fit) || length(fit)!=1)
  stop("'fit' must be a character string")
  pelname<-paste("pel",fit,sep="")
  quaname<-paste("qua",fit,sep="")
  pf<-parent.frame()
  if (!exists(pelname,mode="function",envir=pf)) stop('function "',pelname,'" not found')
  if (!exists(quaname,mode="function",envir=pf)) stop('function "',quaname,'" not found')
  pelfit<-get(pelname,mode="function",envir=pf)
  quafit<-get(quaname,mode="function",envir=pf)

  # Check quantiles

  if (any(f<0 | f>1)) stop("probabilities in 'f' must be between 0 and 1")
  nq<-length(f)

  # Check bound values

  if (any(boundprob<0 | boundprob>1)) stop("probabilities in 'boundprob' must be between 0 and 1")

  # qflocal() is a function that generates a function call to its first
  # argument 'f' with the other arguments of 'qflocal' passed to 'f'.
  # Thus mapply(qflocal,flist,alist) generates a set of calls
  # flist[[i]](alist[[i]]).

  qflocal <- function(f,...) f(...)

  # Simulation loop

  re<-replicate(nrep, {

    # Generate uniform random variates at each site

    if (nocorr) {

      # If there is no inter-site correlation, generate independent
      # uniform samples at each site

      ulist<-lapply(nrec,runif)

    }  else {

      # If there is correlation, generate correlated normal samples
      # each of length 'nmax' ...

      zmat<-matrix(rnorm(nsites*nmax),nsites,nmax)
      zmat<-crossprod(cholcor,zmat)

      # ... transform to uniform ...

      zmat<-pnorm(zmat)

      # ... and retain only the required number of data values at each site

      ulist<-lapply(1:nsites, function(j) zmat[j,1:nrec[j] ])
    }

    # Transform the uniform variates to the required parent distribution.
    # Randomly permute the quantile functions (and their parameters),
    # and, at site i, feed the uniform sample ulist[[i]] through
    # quantile function qfunc[[j]] with parameters para[[j]]
    # where j is the permuted version of i.

    perm<-sample(nsites)
    datlist<-mapply(qflocal, qfunc[perm], ulist, para[perm], SIMPLIFY=FALSE)

    # Rescale the data if necessary, so that all sites have population mean 1,
    # i.e. we permute only the at-site growth curves

    if (!index1) datlist<-mapply("/", datlist, index[perm], SIMPLIFY=FALSE)

    # Rest of loop is just a faster way of computing
    #   xmom <- regsamlmu(datlist)
    #   rgc <- regquant(f,regfit(xmom,fit))
    #   c(xmom[[3]], perm, rgc)

    # Compute L-moments for each site

    xmom<-sapply(datlist,.samlmu,nmom=5)
    sitemeans<-xmom[1,]

    # Compute regional L-moments

    xmom[2,]<-xmom[2,]/xmom[1,]
    xmom[1,]<-1
    rmom<-apply(xmom,1,weighted.mean,w=nrec)

    # Fit the distribution

    rpara<-pelfit(rmom)

    # Regional growth curve estimate at specified quantiles

    rgc<-quafit(f,rpara)

    c(sitemeans,perm,rgc)

  })

  sim.sitemeans<-re[1:nsites,,drop=FALSE]
  sim.perm<-re[nsites+(1:nsites),,drop=FALSE]
  sim.rgc<-re[2*nsites+(1:nq),,drop=FALSE]

  # True at-site quantiles

  trueQ<-mapply(qflocal,qfunc,list(f),para,SIMPLIFY=FALSE)
  trueQ<-matrix(unlist(trueQ),ncol=nsites)

  # True at-site growth curves

  true.asgc<-trueQ/matrix(index,nrow(trueQ),ncol(trueQ),byrow=TRUE)

  # Modified quantile(): if any x value is missing (NA or NaN), set all
  # quantiles equal to NaN.  Avoids error when computing estimation accuracy
  # for an infinite quantile (e.g. in regsimq() when user specifies f=0 and
  # qfunc has no lower bound).
  my.quantile<-function(x,probs,...)
    if (any(is.na(x))) rep(NaN,length(probs)) else quantile(x,probs,type=6)

  # Compute relative RMSE and quantiles of estimated rgc from
  # simulations as an estimator of the at-site growth curve

  sa<-sapply(seq(along=f), function(iq) { # For each quantile F:
    ou<-outer(sim.rgc[iq,,drop=FALSE],true.asgc[iq,,drop=FALSE],"/")
                                         # - matrix of ratios qhat^{[m]}(F)/q_i(F) (F fixed, i varying)
    rr<-sqrt(mean((ou-1)^2))             # - rel. RMSE of qhat as estimator of q_i
    qq<-my.quantile(ou,probs=boundprob)  # - quantiles of the ratio qhat/q_i
    c(rr,qq)
  })
  rel.RMSE<-sa[1,]
  rel.bounds<-t(sa[-1,])
  dimnames(rel.bounds)[[2]]<-boundprob
  relbounds.rgc<-data.frame(f=f,rel.RMSE=rel.RMSE,rel.bound=rel.bounds)

  # Compute relative RMSE and quantiles of quantile estimates
  # for each site

  by.site<-lapply(1:nsites, function(isite) {
    true.asgc.permed<-true.asgc[,sim.perm[isite,] ] # Column j is the growth curve that was used for site 'isite' at repetition j
    rgcratio<-sim.rgc/true.asgc.permed              # Matrix of ratios qhat^{[m]}(F)/q_i(F) (F varying, i fixed)
    meanratio<-sim.sitemeans[isite,]                # Ratio of sample to population mean (sample was generated from distribution with mean 1, so no need to divide by index[isite])
    ratio<-rgcratio*matrix(meanratio,nq,nrep,byrow=TRUE) # Matrix of ratios Qhat_i^{[m]}(F)/q_i(F) (F varying, i fixed)
    rel.RMSE<-sqrt(rowMeans((ratio-1)^2))
    rel.bounds<-t(apply(ratio,1,my.quantile,probs=boundprob))
    dimnames(rel.bounds)[[2]]<-boundprob
    data.frame(f=f,rel.RMSE=rel.RMSE,rel.bound=rel.bounds)
  })

  # sim.rgcratio is a matrix each of whose rows contains, for a single F,
  # 'nrep' realizations of qhat(F)/q_i(F) for varying i
  sim.rgcratio<-NULL
  if (save) {
    sim.rgcratio<-sim.rgc/true.asgc[,rep(1:nsites,length=nrep)]
    rownames(sim.rgcratio)<-f
  }

  out<-list(
    f=f,
    boundprob=boundprob,
    relbounds.rgc=relbounds.rgc,
    relbounds.by.site=by.site,
    sim.rgcratio=sim.rgcratio)

  class(out)<-"regsimq"

  return(out)
}

print.regsimq<-function(x, ...) {
  # Print an object of class "regsimq"
  # Print just the bounds for the regional growth curve, to 3 decimal places
  cat("Regional simulation:", length(x$relbounds.by.site),"sites,",
    ncol(x$sim.rgcratio),"simulations\n")
  cat("Relative RMSE and error bounds for ratio of\nestimated regional growth curve to true at-site growth curve\n")
  print(round(x$relbounds.rgc, 3))
}

regquantbounds<-function(relbounds, rfd) {
## Error bounds for regional growth curve
  if (!inherits(relbounds,"regsimq")) stop("'relbounds' must must be an object of class \"regsimq\"")
  if (!inherits(rfd,"rfd")) stop("'rfd' must be an object of class \"rfd\"")
  if (length(rfd$index)!=length(relbounds$relbounds.by.site))
    warning("regions in 'relbounds' and 'rfd' have different numbers of sites (",
      length(relbounds$relbounds.by.site),",",length(rfd$index),")")
  rgc<-regquant(relbounds$f,rfd)
  RMSE<-abs(rgc)*relbounds$relbounds.rgc$rel.RMSE
  RMSE[rgc==0]<-NaN
  num<-matrix(rgc,length(rgc),length(relbounds$relbounds.rgc)-2,byrow=FALSE)
  denom<-as.matrix(rev(relbounds$relbounds.rgc[-(1:2)]))
  bound<-num/denom
  bound[num>0 & denom<0]<-Inf
  bound[num<=0]<-NA
  out<-cbind(
    f=relbounds$f,
    qhat=rgc,
    RMSE=RMSE,
    bound=as.data.frame(bound))
  boundprob<-rev(1-relbounds$boundprob)
  colnames(out)[-(1:3)]<-paste("bound",boundprob,sep=".")
  rownames(out)<-NULL
  attr(out,"boundprob")<-boundprob
  class(out)<-c("rfdbounds",class(out))
  out
}

sitequantbounds<-function(relbounds, rfd, sitenames, index, seindex, drop=TRUE) {
## Error bounds for quantiles at individual sites
  if (!inherits(relbounds,"regsimq")) stop("'relbounds' must must be an object of class \"regsimq\"")
  if (!inherits(rfd,"rfd")) stop("'rfd' must be an object of class \"rfd\"")
  rgc<-regquant(relbounds$f,rfd)
  boundprob<-rev(1-relbounds$boundprob)
  if (missing(index)) {
    index <- if (missing(sitenames)) rfd$index else rfd$index[sitenames]
    if (any(is.na(index)))
      stop("unable to match 'sitenames' with the sitenames of the region in 'rfd'")
    if (!missing(seindex)) warning("'seindex' ignored when 'index' is missing")
    sitenumbers <- if (missing(sitenames)) seq_along(rfd$index)
      else if (is.numeric(sitenames)) sitenames
      else match(sitenames, names(rfd$index))
    out<-lapply(seq_along(sitenumbers), function(j) {
      isite<-sitenumbers[j]
      relbounds.isite<-relbounds$relbounds.by.site[[isite]]
      Qhat<-index[j]*rgc
      RMSE<-abs(Qhat)*relbounds.isite$rel.RMSE
      RMSE[Qhat==0]<-NaN
      num<-matrix(Qhat,length(Qhat),length(relbounds.isite)-2,byrow=FALSE)
      denom<-as.matrix(rev(relbounds.isite[-(1:2)]))
      bound<-num/denom
      bound[num>0 & denom<0]<-Inf
      bound[num<=0]<-NA
      out.isite<-cbind(
        f=relbounds$f,
        Qhat=Qhat,
        RMSE=RMSE,
        bound=as.data.frame(bound))
      colnames(out.isite)[-(1:3)]<-paste("bound",boundprob,sep=".")
      rownames(out.isite)<-NULL
      attr(out.isite,"boundprob")<-boundprob
      class(out.isite)<-c("rfdbounds",class(out.isite))
      out.isite
    })
    names(out)<-names(rfd$index)[sitenumbers]
  } else {   # 'index' is present
    if (is.null(relbounds$sim.rgcratio))
      stop("'index' cannot be present if 'relbounds$sim.rgcratio' is NULL")
    if (!all(index>=0)) stop("vales of 'index' must not be negative")
    seindex  # Generates a standard error message if 'seindex' is missing
    if (length(seindex)!=length(index))
      stop("'index' and 'seindex' have different lengths")
    if (!all(seindex>=0)) stop("vales of 'seindex' must not be negative")
    if (!missing(sitenames)) {
      if (length(index)!=length(sitenames))
        stop("'index' and 'sitenames' have different lengths")
    }
    nrep<-ncol(relbounds$sim.rgcratio)
    nq  <-nrow(relbounds$sim.rgcratio)
#
    # Modified quantile(): if any x value is missing (NA or NaN), set all quantiles equal to NaN.
    # Avoids error when computing estimation accuracy for an infinite quantile
    # (e.g. in regsimq() when user specifies f=0 and qfunc has no lower bound).
    my.quantile<-function(x,probs,...)
      if (any(is.na(x))) rep(NaN,length(probs)) else quantile(x,probs,type=6)
#
    out<-lapply(seq_along(index),
      function(isite) {
        Qhat<-index[isite]*rgc
        # meanratio: random sample from the distribution of mu_i[hat]/mu_i for this site
        meanratio<-exp(rnorm(nrep)*seindex[isite]/index[isite])
  #     meanratio<-rnorm(nrep,mean=1,sd=seindex[isite]/index[isite]) # alternative
        # Qratio: each row contains a random sample from the distribution of
        # Qhat_i(F)/Q_i(F) for this site
        Qratio<-relbounds$sim.rgcratio*rep(meanratio,each=nq)
        #
        rel.RMSE<-sqrt(rowMeans((Qratio-1)^2))
        RMSE=abs(Qhat)*rel.RMSE
        RMSE[Qhat==0]<-NaN
        #
        rel.bounds<-t(apply(1/Qratio,1,my.quantile,probs=boundprob))
        mult<-matrix(Qhat,nq,ncol(rel.bounds),byrow=FALSE)
        bound<-mult*rel.bounds
        bound[mult>0 & rel.bounds<0]<-Inf
        bound[mult<=0]<-NA
        #
        dimnames(bound)[[2]]<-boundprob
        out.isite<-data.frame(row.names=NULL,
          f=relbounds$f,
          Qhat=Qhat,
          RMSE=RMSE,
          bound=bound)
        attr(out.isite,"boundprob")<-boundprob
        class(out.isite)<-c("rfdbounds",class(out.isite))
        out.isite
      })
    if (!missing(sitenames)) names(out)<-sitenames
    else if (!is.null(names(index))) names(out)<-names(index)
  }
  if (drop && length(out)==1) out<-out[[1]]
  return(out)
}

evplot.rfd<-function(y, ybounds, npoints=101, add=FALSE,
  plim, xlim=c(-2,5), ylim,
  xlab=expression("Reduced variate,  " * -log(-log(italic(F)))),
  ylab="Quantile", rp.axis=TRUE, type="l", lty=c(1,2), col=c(1,1), ...) {
## evplot() method for an object of class "rfd"
## Plots a regional frequency distribution, optionally with error bounds
  if (missing(ybounds)) ybounds<-NULL
  else if (!inherits(ybounds,"rfdbounds")) {
    warning("'ybounds' is not an object of class \"rfdbounds\" -- error bounds not plotted")
    ybounds<-NULL
  }
  if (is.null(ybounds)) {
    if (add)
      evdistq(y$qfunc, npoints=npoints, type=type, lty=lty[1], col=col[1], ...)
    else
      evplot(, qfunc=y$qfunc, npoints=npoints, plim=plim, xlim=xlim, ylim=ylim,
        type=type, xlab=xlab, ylab=ylab, rp.axis=rp.axis,
        lty=lty[1], col=col[1], ...)
    return(invisible())
  }
  #
  mult<-ybounds[[2]]/y$qfunc(ybounds[[1]])
  if (!all.equal(max(mult),min(mult)))
    warning("'ybounds' appears not to be derived from the same distribution as 'y'")
  mult<-median(mult)
  my.qfunc<-function(f) mult*y$qfunc(f)
  if (missing(xlim) && missing(plim)) plim<-range(ybounds[[1]])
  if (missing(ylim)) {
    dat <- if (missing(plim)) my.qfunc(exp(-exp(-xlim))) else my.qfunc(plim)
    dat<-c(dat,0,unlist(ybounds[-(1:3)]))
    ylim<-range(dat[is.finite(dat)])
  }
  #
  if (add)
    evdistq(y$qfunc, npoints=npoints, type=type, lty=lty[1], col=col[1], ...)
  else
    evplot(, qfunc=my.qfunc, npoints=npoints, plim=plim, xlim=xlim, ylim=ylim,
      type=type, xlab=xlab, ylab=ylab, rp.axis=rp.axis,
      lty=lty[1], col=col[1], ...)
  matlines(-log(-log(ybounds[[1]])),ybounds[-(1:3)], type=type,
    lty = if (length(lty)>1) lty[-1] else lty,
    col = if (length(col)>1) col[-1] else col)  # ',...)' ?
  return(invisible())
}
