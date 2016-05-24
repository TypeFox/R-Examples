
###
### main function
###
ecoML <- function(formula, data = parent.frame(), N=NULL, supplement = NULL, 
                  theta.start = c(0,0,1,1,0), fix.rho = FALSE,
                  context = FALSE, sem = TRUE, epsilon=10^(-6),
                  maxit = 1000, loglik = TRUE, hyptest=FALSE, verbose= FALSE) { 

  
  ## getting X and Y
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data=data))
  
  #n.var: total number of parameters involved in the estimation
  #n.par: number of nonstatic paramters need to estimate through EM 
  #       also need SEM
  #ndim: dimension of the multivariate normal distribution

  ndim<-2
  if (context) ndim<-3

  n.var<-2*ndim+ ndim*(ndim-1)/2
 
  n.par<-n.S<-n.var 
  if (context) {
      n.par<-n.var-2
   }

  r12<-NULL
  if (fix.rho) 
     r12<-theta.start[n.par]

  if (!context & fix.rho) n.par<-n.par-1

  flag<-as.integer(context)+2*as.integer(fix.rho)+2^2*as.integer(sem)

  ##checking data
  tmp <- checkdata(X, Y, supplement, ndim)
  bdd <- ecoBD(formula=formula, data=data)
  W1min <- bdd$Wmin[order(tmp$order.old)[1:nrow(tmp$d)],1,1]
  W1max <- bdd$Wmax[order(tmp$order.old)[1:nrow(tmp$d)],1,1]


  n <- tmp$n.samp+tmp$survey.samp+tmp$samp.X1+tmp$samp.X0
  wcol<-ndim
  if (context) {
    wcol<-wcol-1
  }
  inSample.length <- wcol*tmp$n.samp

  #if NCAR and the user did not provide a theta.start
  if (context && (length(theta.start)==5) ) 
    theta.start<-c(0,0,1,1,0,0,0)

  ## Fitting the model via EM  
  res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
            as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
            as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
            as.double(tmp$survey.data),
            as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
            as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
            as.double(W1min), as.double(W1max),
            as.integer(flag),as.integer(verbose),as.integer(loglik),as.integer(hyptest),
            optTheta=rep(-1.1,n.var), pdTheta=double(n.var),
            S=double(n.S+1),inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
            itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
            PACKAGE="eco")

  ##record results from EM
  theta.em<-res$pdTheta
  theta.fisher<-param.trans(theta.em, transformation="Fisher")
  iters.em<-res$itersUsed
  mu.log.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  sigma.log.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  loglike.log.em <- as.double(rep(NA,iters.em))
  nrho<-length(theta.em)-2*ndim
  rho.fisher.em <- matrix(rep(NA,iters.em*nrho),ncol=nrho)
  for(i in 1:iters.em) {
    mu.log.em[i,1:ndim]=res$history[(i-1)*(n.var+1)+(1:ndim)]
    sigma.log.em[i,1:ndim]=res$history[(i-1)*(n.var+1)+ndim+(1:ndim)]
     if (nrho!=0)
    rho.fisher.em[i, 1:nrho]=res$history[(i-1)*(n.var+1)+2*ndim+(1:nrho)]
    loglike.log.em[i]=res$history[(i-1)*(n.var+1)+2*ndim+nrho+1]
  }

  ## In sample prediction of W
  W <- matrix(rep(NA,inSample.length),ncol=wcol)
  for (i in 1:tmp$n.samp)
    for (j in 1:wcol)
      W[i,j]=res$inSample[(i-1)*2+j]

  ## SEM step
  iters.sem<-0

   suff.stat<-res$S
  if (context)
      {
     suff.stat<-rep(0,(n.var+1))
         suff.stat[1]<-mean(logit(c(X,supplement[,3])))
         suff.stat[2:3]<-res$S[1:2]
         suff.stat[4]<-mean((logit(c(X, supplement[,3])))^2)
         suff.stat[5:6]<-res$S[3:4]
         suff.stat[7:8]<-res$S[6:7]
         suff.stat[9]<-res$S[5]
         suff.stat[10]<-res$S[8]
      }


  if (sem) 
  {

    DM <- matrix(rep(NA,n.par*n.par),ncol=n.par)

    res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
              as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
              as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
              as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
              as.integer(flag),as.integer(verbose),as.integer(loglik),as.integer(hyptest),
              res$pdTheta, pdTheta=double(n.var), S=double(n.S+1),
              inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
              itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
              PACKAGE="eco")     
  
    iters.sem<-res$itersUsed
    for(i in 1:n.par)
      for(j in 1:n.par)
        DM[i,j]=res$DMmatrix[(i-1)*n.par+j]


} 

 
   
  if (!context) names(theta.em)<-c("u1","u2","s1","s2","r12")
  if (context) names(theta.em)<-c("ux","u1","u2","sx","s1","s2","r1x","r2x","r12")



  ## output
  res.out<-list(call = mf, Y = Y, X = X, N = N, 
                fix.rho = fix.rho, context = context, sem=sem, epsilon=epsilon,
        theta.em=theta.em, r12=r12, 
               sigma.log = theta.fisher[(ndim+1):(2*ndim)], suff.stat = suff.stat[1:n.S],
                loglik = res$S[n.S+1], iters.em = iters.em, 
                iters.sem = iters.sem, mu.log.em = mu.log.em, 
                sigma.log.em = sigma.log.em,
                rho.fisher.em = rho.fisher.em, loglike.log.em = loglike.log.em,
                W = W)
  
  if (sem) {
    res.out$DM<-DM
#print(dim(data))
# n<-dim(data)[1]
if (!is.null(supplement)) n<-n+dim(supplement)[1]
#cat("n2=", n,"\n")

 res.info<- ecoINFO(theta.em=res.out$theta.em, suff.stat=res.out$suff.stat, DM=res.out$DM, context=context, fix.rho=fix.rho, sem=sem, r12=res.out$r12, n=n)

    res.out$DM<-res.info$DM
    res.out$Icom<-res.info$Icom
    res.out$Iobs<-res.info$Iobs
    res.out$Fmis<-res.info$Fmis
    res.out$Vobs.original<-res.info$Vobs.original
    res.out$Vobs<-res.info$Vobs
    res.out$Iobs<-res.info$Iobs
    res.out$VFmis<-res.info$VFmis
    res.out$Icom.trans<-res.info$Icom.trans
    res.out$Iobs.trans<-res.info$Iobs.trans
    res.out$Fmis.trans<-res.info$Fmis.trans
    res.out$Imiss<-res.info$Imiss
    res.out$Ieigen<-res.info$Ieigen

 res.out$Iobs<-res.info$Iobs
}

  class(res.out) <- "ecoML"
  return(res.out)
}
