

#distance coomputation for distance matrices
gromovl1.intern<-function(d1,d2,n,...){
   #upper bound for rows from lemma 7 in Liebscher(2015a)
   MM<-2*max(abs(d1-d2))
 mmm<-cbind(rep(1:n,n),rep(1:n,each=n))
               mmm<-mmm[which(mmm[,1]<mmm[,2]),]
               ml<-length(mmm[,1])
               rh<-sapply(1:ml,function(p)abs(d1[mmm[p,1],mmm[p,2]]-d2[mmm[p,1],mmm[p,2]]))

  lpd <- initProbGLPK()
  setProbNameGLPK(lpd, "Gromov l1 API")
  setObjDirGLPK(lpd, GLP_MIN)
  addColsGLPK(lpd, n)
  setObjCoefsGLPK(lpd, 1:n,rep(1,n) )
  setColsBndsGLPK(lpd,1:n,lb=rep(0,n),ub=rep(MM,n))
  addRowsGLPK(lpd, ml)
  loadMatrixGLPK(lpd, 2*ml, rep(1:ml,2), c(mmm[,1],mmm[,2]), rep(1,2*ml))
  setRowsBndsGLPK(lpd,1:ml,lb=rh,ub=0*rh+MM)
  setSimplexParmGLPK(102, 3)
  setSimplexParmGLPK(101, 0)
  solveSimplexGLPK(lpd)
  val<-getObjValGLPK(lpd)
  delProbGLPK(lpd)
  return(val)
}


gromovl2.intern<-function(d1,d2,n,...){
  #upper bound for rows from lemma 7 in Liebscher(2015a)
  mmm<-cbind(rep(1:n,n),rep(1:n,each=n))
  mmm<-mmm[which(mmm[,1]<mmm[,2]),]
  ml<-length(mmm[,1])
  cm<-matrix(2,ncol=ml+n,nrow=3)
  cm[1,1:n+ml]<-1
  cm[2,1:ml]<-mmm[,1]
  cm[3,1:ml]<-mmm[,2]
  cm[2,1:n+ml]<-1:n
  rh<-c(sapply(1:ml,function(p)abs(d1[mmm[p,1],mmm[p,2]]-d2[mmm[p,1],mmm[p,2]])),rep(0,n))
  
  return(sqrt((solve.QP.compact(Dmat=2*diag(n),dvec=rep(0,n),Amat=matrix(1,nrow=2,ncol=ml+n),
                                Aind=cm,bvec=rh))$value))
}

gromovlinfinity.intern<-function(d1,d2,n,...){max(abs(d1-d2))/2}

gromovlp.intern<-function(d1,d2,n,p=NULL,...){
  if (missing(p)|is.null(p)) p<-2 
  mmm<-cbind(rep(1:n,n),rep(1:n,each=n))
  mmm<-mmm[which(mmm[,1]<mmm[,2]),]
  ml<-length(mmm[,1])
  rh<-c(sapply(1:ml,function(q)abs(d1[mmm[q,1],mmm[q,2]]-d2[mmm[q,1],mmm[q,2]])),rep(0,n))
  y<-rep(0,n)
  lh<-rbind(t(sapply(1:ml,FUN=function(i){y[mmm[i,1]]<-1;y[mmm[i,2]]<-1;return(y)})),diag(n))
  ofunc<-function(xab,exponent)return(sum(xab^exponent))
  gfunc<-function(xab,exponent) return(exponent*(xab^(exponent-1)))
   
res<-  constrOptim(theta=rep(max(abs(d1-d2)),n),
              f=ofunc,exponent=p,grad=gfunc,
              ci=rh,ui=lh)$value
return(res^(1/p))
}


