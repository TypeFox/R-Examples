#
# Generalised log-linear modelling via EM and Fisher scoring
#
# Library setup
#
packageStartupMessage("This is gllm 0.34")
#
# EM IPF algorithm of Haber AS207
#
emgllm <- function(y,s,X,maxit=1000,tol=0.00001) {
  if (typeof(X)=="language") {
    X<-model.matrix(X)
  }
  X<-cbind(X,double(nrow(X)))
  em<-emgllmfitter(y,s,X,maxit,tol) 
  deviance=2 * sum(em$y * log(em$y / em$f), na.rm=TRUE)
  observed.values<-em$y
  fitted.values<-em$f
  full.table<-em$e
  list(deviance=deviance,
       observed.values=observed.values,
       fitted.values=fitted.values,
       full.table=full.table)
}
emgllmfitter <- function(y,s,X,maxit,tol)
        .C("gllm",
            y = as.double(y),
            ji = as.integer(s-1), # convert indices to C-style, 0:(n-1)
            c = X,
            istop = as.integer(maxit),
            conv = as.double(tol),
            e = double(nrow(X)) + 1, # default inizialization to 1s
            ni = as.integer(nrow(X)),
            nj = as.integer(length(y)),
            nk = as.integer(ncol(X)-1),
            f = double(length(y)),
            PACKAGE="gllm")
#
# Convert scatter vector used by AS207 to matrix for scoring method
#
scatter<- function(y,s) {
  S<-matrix(rep(0,length(y)*length(s)),nrow=length(y))
  for(i in 1:length(s)) if (s[i]<=length(y)) S[s[i],i]<-1
  t(S)
}
#
# Espeland's scoring algorithm
# y=observed table 
# s=scatter vector 
# X=unaugmented design matrix
# m=starting values for complete table
#
scoregllm<-function(y,s,X,m,tol=1e-5) {
#
# Initialize
#
  eps<-0.00001
  call <- match.call()
  formula<-NULL
  if (typeof(X)=="language") {
    formula<-X   
    X<-model.matrix(X)
  }
  X<-t(X)
  S<-scatter(y,s)
  z<-as.vector(S %*% solve(t(S) %*% S, tol=1e-10) %*% y)
  iter<- 0
  olddev<- -1
  deviance<- 0
  b<-qr.solve(t(X),log(m),tol=1e-10)
#
# Main loop
#
  while(abs(olddev-deviance)>tol) {
    iter<-iter+1
    olddev<-deviance
    if (iter>1) m<- as.vector(exp(t(X) %*% b))
    P<- S %*% solve(t(S) %*% diag(m) %*% S, tol=1e-10) %*% t(S) %*% diag(m) 
    A<- P %*% t(X)
    V<- t(A) %*% diag(m) %*% A
    V<- qr.solve(V,tol=1e-10)
    b<- b - V %*% t(A) %*% (m - z)
    f<- as.vector(t(S) %*% m)
    use<-(y>eps & f>eps)
    deviance<- 2.0*sum(y[use]*log(y[use]/f[use]))
  }
  observed.values<-y
  fitted.values<-f
  residuals<- f
  residuals[f>0]<-(y[f>0]-f[f>0])/sqrt(f[f>0])
  full.table<-m
  coefficients<-as.vector(b)
  names(coefficients)<-rownames(X)
  bl<-(rownames(X)=="")
  names(coefficients)[bl]<-paste("beta",1:nrow(X),sep="")[bl]
  se<-sqrt(diag(V))
  df<-length(y)-qr(X)$rank
  res<-list(call=call,formula=formula,
            iter=iter,deviance=deviance,df=df,
            coefficients=coefficients,se=se,V=V,
            observed.values=observed.values,
            fitted.values=fitted.values,
            residuals=residuals,
            full.table=full.table)
  class(res) <- "gllm"
  res
}
#
# Global front end to call either/both routines
#
gllm <- function(y,s,X,method="hybrid",em.maxit=1,tol=0.00001) {
  if (method=="hybrid" || method=="scoring") {
    scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=em.maxit,tol=tol)$full.table)) 
  }else{
    scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=10000,tol=tol)$full.table)) 
  }
}
#
# Summary of results from gllm
#
summary.gllm <- function(object, ...) {
  tab.coef<-data.frame(object$coefficients, 
                       object$se, 
                       exp(object$coefficients),
                       exp(object$coefficients-1.96*object$se),
                       exp(object$coefficients+1.96*object$se),
                       row.names=names(object$coefficients))
  colnames(tab.coef)<-c("Estimate","S.E.","exp(Estimate)",  
                        "Lower 95% CL","Upper 95% CL")
  tab.fitted<-data.frame(object$observed.values, 
                         object$fitted.values, 
                         object$residuals)
  colnames(tab.fitted)<-c("Observed Count","Predicted","Residual")  
  summary<- list()
  summary$call<-object$call
  summary$nobs<-length(object$observed.values)
  summary$nfull<-length(object$full.table)
  summary$mean.cell<-mean(object$observed.values)
  summary$deviance<-object$deviance 
  summary$model.df<-object$df
  summary$coefficients<-tab.coef
  summary$residuals<-tab.fitted
  class(summary) <- "summary.gllm"
  summary
}
#
# Print summary
#
print.summary.gllm <- function(x, digits=NULL, show.residuals = FALSE, ...) {
  if (is.null(digits))
    digits <- options()$digits
  else options(digits=digits)

  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
          "\n\n", sep = "")
  cat("\nNo. cells in observed table: ", x$nobs, "\n",sep="")
  cat("No. cells in complete table: ", x$nfull, "\n",sep="")
  cat("    Mean observed cell size: ", x$mean.cell, "\n",sep="")
  cat("        Model Deviance (df): ", formatC(x$deviance,digits=2,format="f"),
                                        " (", x$model.df, ")\n\n",sep="")

  print(x$coefficients, digits=digits)
  if (show.residuals) {
    cat("\n")
    print(x$residuals, digits=digits)
  }
}
#
# Bootstrap a contingency table.  Sampling zeroes augmented by 1/2.
#
boot.table <- function(y,strata=NULL) {  
  ynew<-rep(0.5, length(y))
  if (is.null(strata)) {
    tab.ynew<-table(sample(rep(1:length(y),y),replace=TRUE))
    ynew[as.integer(names(tab.ynew))]<-tab.ynew
  }else{
    s<-as.integer(strata)
    for(i in unique(s)) {
      idx<-s==i
      tab.ynew<-table(sample(rep((1:length(y))[idx],y[idx]),replace=TRUE))
      ynew[as.integer(names(tab.ynew))]<-tab.ynew
    }
  }
  ynew
}
#
# Bootstrap a GLLM returning the full fitted tables
#
boot.gllm <- function(y,s,X,method="hybrid",em.maxit=1,tol=0.00001,
                      strata=NULL,R=200) {
  if (method=="hybrid" || method=="scoring") {
    f0<-scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=em.maxit,tol=tol)$full.table))$full.table
  }else{
    f0<-scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=10000,tol=tol)$full.table))$full.table
  }
  result<-as.matrix(t(f0))
  cat("It",0,":",y,"\n")
  for(i in 1:R) {
    ynew<-boot.table(y,strata=strata)
    cat("It",i,":",ynew,"\n")
    result<-rbind(result, t(scoregllm(ynew,s,X,as.array(emgllm(ynew,s,X,
                            maxit=em.maxit,tol=tol)$full.table))$full.table))
  }
  result
}
#
# After anova.multinom
#
anova.gllm <- function(object, ..., test = c("Chisq", "none"))
{
  modelname<-unlist(strsplit(deparse(match.call()),"[(),]")) 
  modelname<-gsub("[() ]","",modelname)
  modelname<-gsub("object=","",modelname)
  modelname<-modelname[-c(1,grep("=",modelname))]

  test <- match.arg(test)
  dots <- list(...)

  mlist <- list(object, ...)
  nt <- length(mlist)
  dflis <- sapply(mlist, function(x) x$df)
  s <- order(dflis, decreasing=TRUE)
  mlist <- mlist[s]
  if(any(!sapply(mlist, inherits, "gllm"))) {
    stop("not all objects are of class `gllm'")
  }
  mds <- modelname[s]
  dfs <- dflis[s]
  lls <- sapply(mlist, function(x) x$deviance)
  tss <- ""
  if (nt>1) tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
  out <- data.frame(Model = mds, Resid.df = dfs,
                    Deviance = lls, Pr.Fit = 1-pchisq(lls,dfs), 
                    Test = tss, Df = df, LRtest = x2,
                    Prob = pr)
  names(out) <- c("Model", "Resid. df", "Resid. Dev", "Pr(GOFChi)",
                  "Test", "   Df", "LR stat.", "Pr(Chi)")
  if(test=="none") out <- out[, 1:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of Loglinear Models\n")
  out
}
#
# Test for allelic association following Aston and Wilson.
# 2 locus version
#
ld2 <- function(locus1, locus2=NULL) {
  all.possible.genos <- function(genotypes,gtp.sep="/") {
    alleles<- sort(unique(unlist(strsplit(genotypes[!is.na(genotypes)], gtp.sep))))
    gtp<-cbind(rep(alleles,seq(length(alleles),1,-1)),
               rev(alleles)[sequence(seq(length(alleles),1,-1))])
    paste(gtp[,1],gtp[,2],sep=gtp.sep)
  }
  genofactor <- function(genotypes, gtp.sep="/") {
    factor(genotypes, levels=all.possible.genos(genotypes))
  }
  ng.to.nall<-function(ngenos) (sqrt(1+8*ngenos)-1)/2

  if (is.null(locus2)) {
    if (is.table(locus1) && length(dim(locus1))==2) {
      ld.table<-locus1
    }else if (is.data.frame(locus1)) {
      ld.table<-table(genofactor(locus1[,1]), genofactor(locus1[,2]))
    }
  }else{
    ld.table<-table(genofactor(locus1), genofactor(locus2))
  }
  siz<-dim(ld.table)
  if (length(siz)==2 && siz[1]>1 && all(siz>1)) {
    nall1<-ng.to.nall(siz[1])
    nall2<-ng.to.nall(siz[2])
    m0<-ld2.model(nall1,nall2,"~a1+a2")
    m1<-ld2.model(nall1,nall2,"~a1+a2+d")
    m2<-ld2.model(nall1,nall2,"~a1+a2+p1")
    m3<-ld2.model(nall1,nall2,"~a1+a2+p1+p2")
    m4<-ld2.model(nall1,nall2,"~a1+a2+p1+p2+d")
    m0<-gllm(c(t(ld.table)),m0$s,m0$X,method="hybrid",em.maxit=1,tol=0.00001) 
    m1<-gllm(c(t(ld.table)),m1$s,m1$X,method="hybrid",em.maxit=1,tol=0.00001) 
    m2<-gllm(c(t(ld.table)),m2$s,m2$X,method="hybrid",em.maxit=1,tol=0.00001) 
    m3<-gllm(c(t(ld.table)),m3$s,m3$X,method="hybrid",em.maxit=1,tol=0.00001) 
    m4<-gllm(c(t(ld.table)),m4$s,m4$X,method="hybrid",em.maxit=1,tol=0.00001) 
    res<-list(m0=m0,m1=m1,m2=m2,m3=m3,m4=m4)
    class(res)<-"ld2"
    res
  }else{
    cat("Must have two polymorphic loci\n")
  }
}
print.ld2 <- function(x,...) {
  cat("\nAssuming HWE\n")
  print(summary(x$m1))
  cat("\n")
  print(anova(x$m0,x$m1))
  cat("\nModelling HWD\n")
  print(summary(x$m4))
  cat("\n")
  print(anova(x$m0,x$m2,x$m3,x$m4))
}
#
# Produce the filter and design matrices.
# Geno is not used directly, but documents the ordering of the 
# observed table of genotype counts
#
ld2.model <- function(nall1, nall2, formula="~a1+a2+p1+p2+d") {
  mkgeno<-function(nall) {
    lab<-t(outer(1:nall, 1:nall, paste,sep="/"))
    lab[lower.tri(lab,diag=TRUE)]
  }
  ld.terms<-function(nall1, nall2, formula=formula) {
    nter<-cumsum(c(1,nall1-1,nall2-1, (nall1-1)*(nall1-1),
                   (nall2-1)*(nall2-1), (nall1-1)*(nall2-1)))
    ter<-pmatch(unlist(strsplit(formula,"[~+]")),c("a1","a2","p1","p2","d"))
    ter<-ter[!is.na(ter)]
    res<-1
    for(i in ter) res<-append(res,seq(nter[i]+1, nter[i+1]))
    res           
  }
  ng1<-nall1*(nall1+1)/2
  ng2<-nall2*(nall2+1)/2
  gam<-nall1*nall1*nall2*nall2
  gen<-ng1*ng2
  mod<-nall1+nall2-1+ 
       (nall1-1)*(nall1-1)+(nall2-1)*(nall2-1)+ 
       (nall1-1)*(nall2-1)
  Geno<-matrix(1:gen, nrow=ng1, byrow=TRUE)
  rownames(Geno)<-mkgeno(nall1)
  colnames(Geno)<-mkgeno(nall2)
  S<-array(0, dim=c(nall1,nall1,nall2,nall2))
  X<-matrix(0, nrow=gam, ncol=mod)
  colnames(X)<-c("Int",paste("a1",2:nall1,sep="."),
                       paste("a2",2:nall2,sep="."),
                       paste("p1",outer(2:nall1,2:nall1,paste,sep=""),sep="."),
                       paste("p2",outer(2:nall2,2:nall2,paste,sep=""),sep="."),
                       paste("d",outer(2:nall1,2:nall2,paste,sep=""),sep="."))
#
# Produce all nall1^2 * nall2^2 gametes and map onto genotypes
#
  k<-0
  for(i in 1:nall1) for(i2 in 1:i)
  for(j in 1:nall2) for(j2 in 1:j) {
    k<-k+1
    S[i,i2,j,j2]<-S[i,i2,j2,j]<-S[i2,i,j,j2]<-S[i2,i,j2,j]<-k
  }
  S<-matrix(S, nrow=nall1*nall1, ncol=nall2*nall2)
  rownames(S)<-outer(1:nall1,1:nall1,paste,sep="/") 
  colnames(S)<-outer(1:nall2,1:nall2,paste,sep="/") 
#
# Then the design matrix
#
  k<-0
  hap<-rep(" ",gam)
  for(i in 1:nall1) for(i2 in 1:nall1)
  for(j in 1:nall2) for(j2 in 1:nall2) {
    k<-k+1
    a1<-rep(0,nall1)
    a1[i]<-a1[i]+1
    a1[i2]<-a1[i2]+1
    p1<-matrix(0, nrow=nall1, ncol=nall1)
    p1[i,i2]<-p1[i,i2]+1
    a2<-rep(0,nall2)
    a2[j]<-a2[j]+1
    a2[j2]<-a2[j2]+1
    p2<-matrix(0, nrow=nall2, ncol=nall2)
    p2[j,j2]<-p2[j,j2]+1
    d<-matrix(0, nrow=nall1, ncol=nall2)
    d[i, j2]<-d[i, j2]+1
    d[i2, j]<-d[i2, j]+1
    X[k,]<-c(1,a1[-1],a2[-1],c(p1[-1,-1]), c(p2[-1,-1]), c(d[-1,-1]))
    hap[k]<-paste(i,j,";",i2,j2,sep="")
  }
  rownames(X)<-hap
  s<-c(t(S))
  names(s)<-hap
  list(Geno=Geno, s=s, X=X[,ld.terms(nall1,nall2,formula)])
}
