##
## Tests for contingency tables
##


svychisq<-function(formula, design,...) UseMethod("svychisq",design)


svychisq.survey.design<-function(formula, design,
                   statistic=c("F","Chisq","Wald","adjWald","lincom","saddlepoint"),
                   na.rm=TRUE,...){
  if (ncol(attr(terms(formula),"factors"))>2)
    stop("Only 2-way tables at the moment")
  statistic<-match.arg(statistic)
  
  ##if(!is.null(design$postStrata))
  ##  warning("Post-stratification not implemented")
  
  rows<-formula[[2]][[2]]
  cols<-formula[[2]][[3]]
  rowvar<-unique(design$variables[,as.character(rows)])
  colvar<-unique(design$variables[,as.character(cols)])
  returnNA<-FALSE
  if ((any(is.na(rowvar),is.na(colvar)))){
      rowvar<-na.omit(rowvar)
      colvar<-na.omit(colvar)
      returnNA<-!na.rm
  }
  nr<-length(rowvar)
  nc<-length(colvar)
  
  fsat<-eval(bquote(~interaction(factor(.(rows)),factor(.(cols)))-1))
  mm<-model.matrix(fsat,model.frame(fsat, design$variables,na.action=na.pass))
  N<-nrow(mm)
  nu <- length(unique(design$cluster[,1]))-length(unique(design$strata[,1]))


  pearson<- suppressWarnings(chisq.test(svytable(formula,design,Ntotal=N),
                       correct=FALSE))

  
  mf1<-expand.grid(rows=1:nr,cols=1:nc)
  X1<-model.matrix(~factor(rows)+factor(cols),mf1)
  X12<-model.matrix(~factor(rows)*factor(cols),mf1)

  
  if(statistic %in% c("Wald", "adjWald")){
    frow<-eval(bquote(~factor(.(rows))-1))
    fcol<-eval(bquote(~factor(.(cols))-1))
    mr<-model.matrix(frow, model.frame(frow,design$variables, na.action=na.pass))
    mc<-model.matrix(fcol, model.frame(fcol,design$variables, na.action=na.pass))
    one<-rep(1,NROW(mc))
    cells<-svytotal(~mm+mr+mc+one,design,na.rm=TRUE)

    Jcb <- cbind(diag(nr*nc),
                 -outer(mf1$rows,1:nr,"==")*rep(cells[(nr*nc)+nr+1:nc]/cells[(nr*nc)+nr+nc+1],each=nr),
                 -outer(mf1$cols,1:nc,"==")*cells[(nr*nc)+1:nr]/cells[(nr*nc)+nr+nc+1],
                 as.vector(outer(cells[(nr*nc)+1:nr],cells[(nr*nc+nr)+1:nc])/cells[(nr*nc)+nr+nc+1]^2))

    Y<-cells[1:(nc*nr)]-as.vector(outer(cells[(nr*nc)+1:nr],cells[(nr*nc+nr)+1:nc]))/cells[(nr*nc)+nr+nc+1]
    V<-Jcb%*%attr(cells,"var")%*%t(Jcb)
    use<-as.vector(matrix(1:(nr*nc),nrow=nr,ncol=nc)[-1,-1])
    waldstat<-Y[use]%*%solve(V[use,use],Y[use])
    if (statistic=="Wald"){
      waldstat<-waldstat/((nc-1)*(nr-1))
      numdf<-(nc-1)*(nr-1)
      denomdf<-nu
    } else {
      numdf<-(nr-1)*(nc-1)
      denomdf<-(nu-numdf+1)
      waldstat <- waldstat*denomdf/(numdf*nu)
    }
    if (returnNA){
      pearson$statistic<-NA
      pearson$parameter<-c(ndf=numdf,ddf=denomdf)
      pearson$p.value<-NA
      attr(pearson$statistic,"names")<-"F"
      pearson$data.name<-deparse(sys.call(-1))
      pearson$method<-"Design-based Wald test of association"
    } else {
      pearson$statistic<-waldstat
      pearson$parameter<-c(ndf=numdf,ddf=denomdf)
      pearson$p.value<-pf(pearson$statistic, numdf, denomdf, lower.tail=FALSE)
      attr(pearson$statistic,"names")<-"F"
      pearson$data.name<-deparse(sys.call(-1))
      pearson$method<-"Design-based Wald test of association"
    } 
    return(pearson)
  }
  
  mean2<-svymean(mm,design,na.rm=TRUE)



  
  Cmat<-qr.resid(qr(X1),X12[,-(1:(nr+nc-1)),drop=FALSE])
  Dmat <- diag(mean2)
  iDmat<- diag(ifelse(mean2==0,0,1/mean2))
  Vsrs <- (Dmat - outer(mean2,mean2))/N
  V <- attr(mean2,"var")
  denom<- t(Cmat) %*% (iDmat/N) %*% Cmat
  numr<-t(Cmat)%*% iDmat %*% V %*% iDmat %*% Cmat
  Delta<-solve(denom,numr)
  d0<- sum(diag(Delta))^2/(sum(diag(Delta%*%Delta)))
  
  warn<-options(warn=-1) ## turn off the small-cell count warning.
  pearson<- chisq.test(svytable(formula,design,Ntotal=N),
                       correct=FALSE)
  options(warn)
  
  if (match.arg(statistic)=="F"){
    pearson$statistic<-pearson$statistic/sum(diag(Delta))
    pearson$p.value<-pf(pearson$statistic, d0, d0*nu, lower.tail=FALSE)
    attr(pearson$statistic,"names")<-"F"
    pearson$parameter<-c(ndf=d0,ddf=d0*nu)
    pearson$method<-"Pearson's X^2: Rao & Scott adjustment"
  }  else if (match.arg(statistic)=="lincom") {
    pearson$p.value<-pFsum(pearson$statistic, rep(1,ncol(Delta)), eigen(Delta,only.values=TRUE)$values,
                              lower.tail=FALSE,method="integration",ddf=d0*nu)
    pearson$parameter<-NULL
    pearson$method<-"Pearson's X^2: asymptotic exact distribution"
  } else if  (match.arg(statistic)=="saddlepoint") {
    pearson$p.value<-pFsum(pearson$statistic, rep(1,ncol(Delta)), eigen(Delta,only.values=TRUE)$values,
                              lower.tail=FALSE,method="saddlepoint",ddf=d0*nu)
    pearson$parameter<-NULL
    pearson$method<-"Pearson's X^2: saddlepoint approximation"
  } else{
    pearson$p.value<-pchisq(pearson$statistic/mean(diag(Delta)),
                            df=NCOL(Delta),lower.tail=FALSE)
    pearson$parameter<-c(df=NCOL(Delta))
    pearson$method<-"Pearson's X^2: Rao & Scott adjustment"
  }

  if (returnNA){
    pearson$statistic<-NA
    pearson$p.value<-NA
  }
  
  pearson$data.name<-deparse(sys.call(-1))

  pearson
  
}

svychisq.twophase<-function(formula, design,
                   statistic=c("F","Chisq","Wald","adjWald","lincom","saddlepoint"),
                   na.rm=TRUE,...){
  if (ncol(attr(terms(formula),"factors"))>2)
    stop("Only 2-way tables at the moment")
  statistic<-match.arg(statistic)
  
  ##if(!is.null(design$postStrata))
  ##  warning("Post-stratification not implemented")
  
  rows<-formula[[2]][[2]]
  cols<-formula[[2]][[3]]
  dat<-design$phase1$sample$variables
  rowvar<-unique(dat[,as.character(rows)])
  colvar<-unique(dat[,as.character(cols)])
  returnNA<-FALSE
  if ((any(is.na(rowvar),is.na(colvar)))){
      rowvar<-na.omit(rowvar)
      colvar<-na.omit(colvar)
      returnNA<-!na.rm
  }
  nr<-length(rowvar)
  nc<-length(colvar)
  
  fsat<-eval(bquote(~interaction(factor(.(rows)),factor(.(cols)))-1))
  mm<-model.matrix(fsat,model.frame(fsat, dat,na.action=na.pass))
  N<-nrow(mm)
  nu <- length(unique(design$phase2$cluster[,1]))-length(unique(design$phase2$strata[,1]))


  pearson<- suppressWarnings(chisq.test(svytable(formula,design,Ntotal=N),
                       correct=FALSE))

  
  mf1<-expand.grid(rows=1:nr,cols=1:nc)
  X1<-model.matrix(~factor(rows)+factor(cols),mf1)
  X12<-model.matrix(~factor(rows)*factor(cols),mf1)

  
  if(statistic %in% c("Wald", "adjWald")){
    frow<-eval(bquote(~factor(.(rows))-1))
    fcol<-eval(bquote(~factor(.(cols))-1))
    mr<-model.matrix(frow, model.frame(frow,dat, na.action=na.pass))
    mc<-model.matrix(fcol, model.frame(fcol,dat, na.action=na.pass))
    one<-rep(1,NROW(mc))
    cells<-svytotal(~mm+mr+mc+one,design,na.rm=TRUE)

    Jcb <- cbind(diag(nr*nc),
                 -outer(mf1$rows,1:nr,"==")*rep(cells[(nr*nc)+nr+1:nc]/cells[(nr*nc)+nr+nc+1],each=nr),
                 -outer(mf1$cols,1:nc,"==")*cells[(nr*nc)+1:nr]/cells[(nr*nc)+nr+nc+1],
                 as.vector(outer(cells[(nr*nc)+1:nr],cells[(nr*nc+nr)+1:nc])/cells[(nr*nc)+nr+nc+1]^2))

    Y<-cells[1:(nc*nr)]-as.vector(outer(cells[(nr*nc)+1:nr],cells[(nr*nc+nr)+1:nc]))/cells[(nr*nc)+nr+nc+1]
    V<-Jcb%*%attr(cells,"var")%*%t(Jcb)
    use<-as.vector(matrix(1:(nr*nc),nrow=nr,ncol=nc)[-1,-1])
    waldstat<-Y[use]%*%solve(V[use,use],Y[use])
    if (statistic=="Wald"){
      waldstat<-waldstat/((nc-1)*(nr-1))
      numdf<-(nc-1)*(nr-1)
      denomdf<-nu
    } else {
      numdf<-(nr-1)*(nc-1)
      denomdf<-(nu-numdf+1)
      waldstat <- waldstat*denomdf/(numdf*nu)
    }
    if (returnNA){
      pearson$statistic<-NA
      pearson$parameter<-c(ndf=numdf,ddf=denomdf)
      pearson$p.value<-NA
      attr(pearson$statistic,"names")<-"F"
      pearson$data.name<-deparse(sys.call(-1))
      pearson$method<-"Design-based Wald test of association"
    } else {
      pearson$statistic<-waldstat
      pearson$parameter<-c(ndf=numdf,ddf=denomdf)
      pearson$p.value<-pf(pearson$statistic, numdf, denomdf, lower.tail=FALSE)
      attr(pearson$statistic,"names")<-"F"
      pearson$data.name<-deparse(sys.call(-1))
      pearson$method<-"Design-based Wald test of association"
    } 
    return(pearson)
  }
  
  mean2<-svymean(mm,design,na.rm=TRUE)



  
  Cmat<-qr.resid(qr(X1),X12[,-(1:(nr+nc-1)),drop=FALSE])
  Dmat <- diag(mean2)
  iDmat<- diag(ifelse(mean2==0,0,1/mean2))
  Vsrs <- (Dmat - outer(mean2,mean2))/N
  V <- attr(mean2,"var")
  denom<- t(Cmat) %*% (iDmat/N) %*% Cmat
  numr<-t(Cmat)%*% iDmat %*% V %*% iDmat %*% Cmat
  Delta<-solve(denom,numr)
  d0<- sum(diag(Delta))^2/(sum(diag(Delta%*%Delta)))
  
  warn<-options(warn=-1) ## turn off the small-cell count warning.
  pearson<- chisq.test(svytable(formula,design,Ntotal=N),
                       correct=FALSE)
  options(warn)
  
  if (match.arg(statistic)=="F"){
    pearson$statistic<-pearson$statistic/sum(diag(Delta))
    pearson$p.value<-pf(pearson$statistic, d0, d0*nu, lower.tail=FALSE)
    attr(pearson$statistic,"names")<-"F"
    pearson$parameter<-c(ndf=d0,ddf=d0*nu)
    pearson$method<-"Pearson's X^2: Rao & Scott adjustment"
  }  else if (match.arg(statistic)=="lincom") {
    pearson$p.value<-pFsum(pearson$statistic, rep(1,ncol(Delta)), eigen(Delta,only.values=TRUE)$values,
                              lower.tail=FALSE,method="integration",ddf=d0*nu)
    pearson$parameter<-NULL
    pearson$method<-"Pearson's X^2: asymptotic exact distribution"
  } else if  (match.arg(statistic)=="saddlepoint") {
    pearson$p.value<-pFsum(pearson$statistic, rep(1,ncol(Delta)), eigen(Delta,only.values=TRUE)$values,
                              lower.tail=FALSE,method="saddlepoint",ddf=d0*nu)
    pearson$parameter<-NULL
    pearson$method<-"Pearson's X^2: saddlepoint approximation"
  } else{
    pearson$p.value<-pchisq(pearson$statistic/mean(diag(Delta)),
                            df=NCOL(Delta),lower.tail=FALSE)
    pearson$parameter<-c(df=NCOL(Delta))
    pearson$method<-"Pearson's X^2: Rao & Scott adjustment"
  }

  if (returnNA){
    pearson$statistic<-NA
    pearson$p.value<-NA
  }
  
  pearson$data.name<-deparse(sys.call(-1))

  pearson
  
}


svychisq.svyrep.design<-function(formula, design,
                   statistic=c("F","Chisq","Wald","adjWald","lincom","saddlepoint"),
                   na.rm=TRUE,...){
  if (ncol(attr(terms(formula),"factors"))>2)
    stop("Only 2-way tables at the moment")
  statistic<-match.arg(statistic)
    
  rows<-formula[[2]][[2]]
  cols<-formula[[2]][[3]]
  rowvar<-unique(design$variables[,as.character(rows)])
  colvar<-unique(design$variables[,as.character(cols)])
  returnNA<-FALSE
  if ((any(is.na(rowvar),is.na(colvar)))){
      rowvar<-na.omit(rowvar)
      colvar<-na.omit(colvar)
      returnNA<-!na.rm
  }
  nr<-length(rowvar)
  nc<-length(colvar)
  
  fsat<-eval(bquote(~interaction(factor(.(rows)),factor(.(cols)))-1))
  mm<-model.matrix(fsat,model.frame(fsat, design$variables,na.action=na.pass))
  N<-nrow(mm)
  nu <- degf(design)


  pearson<- suppressWarnings(chisq.test(svytable(formula,design,Ntotal=N),
                       correct=FALSE))

  
  mf1<-expand.grid(rows=1:nr,cols=1:nc)
  X1<-model.matrix(~factor(rows)+factor(cols),mf1)
  X12<-model.matrix(~factor(rows)*factor(cols),mf1)

  
  if(statistic %in% c("Wald", "adjWald")){
    frow<-eval(bquote(~factor(.(rows))-1))
    fcol<-eval(bquote(~factor(.(cols))-1))
    mr<-model.matrix(frow, model.frame(frow,design$variables, na.action=na.pass))
    mc<-model.matrix(fcol, model.frame(fcol,design$variables, na.action=na.pass))
    one<-rep(1,NROW(mc))
    cells<-svytotal(~mm+mr+mc+one,design,na.rm=TRUE)

    Jcb <- cbind(diag(nr*nc),
                 -outer(mf1$rows,1:nr,"==")*rep(cells[(nr*nc)+nr+1:nc]/cells[(nr*nc)+nr+nc+1],each=nr),
                 -outer(mf1$cols,1:nc,"==")*cells[(nr*nc)+1:nr]/cells[(nr*nc)+nr+nc+1],
                 as.vector(outer(cells[(nr*nc)+1:nr],cells[(nr*nc+nr)+1:nc])/cells[(nr*nc)+nr+nc+1]^2))

    Y<-cells[1:(nc*nr)]-as.vector(outer(cells[(nr*nc)+1:nr],cells[(nr*nc+nr)+1:nc]))/cells[(nr*nc)+nr+nc+1]
    V<-Jcb%*%attr(cells,"var")%*%t(Jcb)
    use<-as.vector(matrix(1:(nr*nc),nrow=nr,ncol=nc)[-1,-1])
    waldstat<-Y[use]%*%solve(V[use,use],Y[use])
    if (statistic=="Wald"){
      waldstat<-waldstat/((nc-1)*(nr-1))
      numdf<-(nc-1)*(nr-1)
      denomdf<-nu
    } else {
      numdf<-(nr-1)*(nc-1)
      denomdf<-(nu-numdf+1)
      waldstat <- waldstat*denomdf/(numdf*nu)
    }
    if (returnNA){
      pearson$statistic<-NA
      pearson$parameter<-c(ndf=numdf,ddf=denomdf)
      pearson$p.value<-NA
      attr(pearson$statistic,"names")<-"F"
      pearson$data.name<-deparse(sys.call(-1))
      pearson$method<-"Design-based Wald test of association"
    } else {
      pearson$statistic<-waldstat
      pearson$parameter<-c(ndf=numdf,ddf=denomdf)
      pearson$p.value<-pf(pearson$statistic, numdf, denomdf, lower.tail=FALSE)
      attr(pearson$statistic,"names")<-"F"
      pearson$data.name<-deparse(sys.call(-1))
      pearson$method<-"Design-based Wald test of association"
    } 
    return(pearson)
  }
  
  mean2<-svymean(mm,design,na.rm=TRUE)



  
  Cmat<-qr.resid(qr(X1),X12[,-(1:(nr+nc-1)),drop=FALSE])
  Dmat <- diag(mean2)
  iDmat<- diag(ifelse(mean2==0,0,1/mean2))
  Vsrs <- (Dmat - outer(mean2,mean2))/N
  V <- attr(mean2,"var")
  denom<- t(Cmat) %*% (iDmat/N) %*% Cmat
  numr<-t(Cmat)%*% iDmat %*% V %*% iDmat %*% Cmat
  Delta<-solve(denom,numr)
  d0<- sum(diag(Delta))^2/(sum(diag(Delta%*%Delta)))
  
  warn<-options(warn=-1) ## turn off the small-cell count warning.
  pearson<- chisq.test(svytable(formula,design,Ntotal=N),
                       correct=FALSE)
  options(warn)
  
  if (match.arg(statistic)=="F"){
    pearson$statistic<-pearson$statistic/sum(diag(Delta))
    pearson$p.value<-pf(pearson$statistic, d0, d0*nu, lower.tail=FALSE)
    attr(pearson$statistic,"names")<-"F"
    pearson$parameter<-c(ndf=d0,ddf=d0*nu)
  }  else if (match.arg(statistic)=="lincom") {
    pearson$p.value<-pchisqsum(pearson$statistic, rep(1,ncol(Delta)), eigen(Delta,only.values=TRUE)$values,
                              lower.tail=FALSE,method="integration")
    pearson$parameter<-NULL
    pearson$method<-"Pearson's X^2: asymptotic exact distribution"
  }else if  (match.arg(statistic)=="saddlepoint") {
    pearson$p.value<-pchisqsum(pearson$statistic, rep(1,ncol(Delta)), eigen(Delta,only.values=TRUE)$values,
                              lower.tail=FALSE,method="saddlepoint")
    pearson$parameter<-NULL
    pearson$method<-"Pearson's X^2: saddlepoint approximation"
  }  else {
    pearson$p.value<-pchisq(pearson$statistic/mean(diag(Delta)),
                               df=NCOL(Delta),lower.tail=FALSE)
    pearson$parameter<-c(df=NCOL(Delta))
  }

  if (returnNA){
    pearson$statistic<-NA
    pearson$p.value<-NA
  }
  
  pearson$data.name<-deparse(sys.call(-1))
  pearson$method<-"Pearson's X^2: Rao & Scott adjustment"
  pearson
  
}



summary.svreptable<-function(object,...){
  object
}

summary.svytable<-function(object, statistic=c("F","Chisq","Wald","adjWald","lincom","saddlepoint"),...){

  statistic<-match.arg(statistic)
  call<-attr(object, "call")
  ff<-call$formula
  
  if (is.null(environment(ff)))
    env<-parent.frame()
  else
    env<-environment(ff)

  ff<-delete.response(ff)
      
  test<-eval(bquote(svychisq(.(ff), design=.(call$design),
                             statistic=.(statistic))), env)

  rval<-list(table=object,statistic=test)
  class(rval)<-"summary.svytable"
  rval
}

print.summary.svytable<-function(x,digits=0,...){
  print(round(x$table,digits))
  print(x$statistic,...)
}
