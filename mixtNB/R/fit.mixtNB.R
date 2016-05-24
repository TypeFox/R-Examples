#################
# EM algorithm  #
#################


fit.mixtNB<-function(y,cr,K,it=200,eps=1e-5,init=NULL,seme=1,filter=TRUE,quiet=FALSE)
{
  ptm <- proc.time()
  set.seed(seme)
  
  p=dim(y)[1]
  gname=(1:p)
  a.max=500
  
  if (filter) {
  out=filter.em(y)
  y=out$y
  n.rm=out$n.rm
  gname=out$gname
  if (!quiet) cat(n.rm,"genes has been filtered because they contains too small number of reads across the experiments.\n")
  
  } else n.rm=0
 
  d<-length(unique(cr))
  p<-dim(y)[1] # number of units
  n<-dim(y)[2] # total number of replicates (number of columns)
  lambda<-matrix(0,p,d) # matrix for the estimation of parameters lambda_ij
  a<-matrix(0,K,1)
  w<-matrix(0,K,1)
  y<-as.matrix(y)
  nj=n/d
  
  #############################
  # Initialization
  #############################


  
  if (is.null(init$w)) {w<-rep(1/K,K)
                        w<-w/sum(w)} else w<-init$w
  for (j in 1:d) lambda[,j]<-rowMeans(y[,cr==j,drop=FALSE])
  if (is.null(init$a))  a<-round(seq(0.1,a.max-5,length.out=K),2) else a<-init$a
  
  likelihood<-NULL
  ratio = 1000
  lik<--10000000
  
  f.y.z<-array(0,c(p,n,K))
  f.y.zProd<-matrix(0,p,K)
  post<-matrix(0,p,K)
  f.z.y<-matrix(0,p,K) # posterior
  a_tilde<-array(0,c(p,n,K))
  b_tilde<-array(0,c(p,n,K))
  E.u.yz<-array(0,c(p,n,K))
  E.logu.yz<-array(0,c(p,n,K))
 
  h=0
  
  #################
  ### EM begins ###
  #################
  
  while ((h < it) & (ratio > eps)) {

    f.z.y_OLD<-f.z.y
    
    #################
    #### E-STEP #####
    #################
    
    for (k in 1:K)  {for (l in 1:n) {
                                  f.y.z[,l,k]<-dnbinom(y[,l],mu=lambda[,cr[l]],size=a[k],log=TRUE)
                                  a_tilde[,l,k]<-y[,l]+a[k]
                                  b_tilde[,l,k]<-lambda[,cr[l]]+a[k]
                                      }
                     f.y.zProd[,k]<-exp(rowSums(f.y.z[,,k,drop=FALSE]))
                     f.y.zProd[,k]<-ifelse(is.na(f.y.zProd[,k]),10^(-300),f.y.zProd[,k,drop=FALSE])
                     post[,k]<-f.y.zProd[,k,drop=FALSE]*w[k]
                      }
    
  post[rowSums(post)==0,]=10^(-300)
                         
    
  f.z.y<-post/matrix(rowSums(post),p,K)
  E.u.yz<-a_tilde/b_tilde
  E.logu.yz<-digamma(a_tilde)-log(b_tilde)

  E.u.yz<-ifelse(is.na(E.u.yz),mean(E.u.yz,na.rm=TRUE), E.u.yz)
  E.logu.yz<-ifelse(is.na(E.logu.yz),mean(E.logu.yz,na.rm=TRUE),E.logu.yz)
      
    #################
    #### M-STEP #####
    #################
   
     for (k in 1:K) {
        a[k]=nlminb(start=a[k],objective=fun.a,lower =0.01, upper =500, n=n,p=p,
                    E.u.yz=E.u.yz[,,k], E.logu.yz=E.logu.yz[,,k,drop=FALSE],f.z.y=f.z.y[,k],control=list(iter.max=5))$par
        
        w[k]<-mean(f.z.y[,k,drop=FALSE])
     }
    
   temp<-sum(log(f.y.zProd%*%w))
   likelihood = c(likelihood, temp)
   ratio = (temp - lik)/abs(lik)
   lik = temp
   
    h<-h+1
  }
  
  if (!quiet) cat("The EM converges after",h,"iterations.\n")
  
  npar<-2*K-1+p*d

  AIC<--2*likelihood[length(likelihood)]+2*npar
  BIC<--2*likelihood[length(likelihood)]+log(p)*npar
  disp=matrix(rowSums(f.z.y/matrix(a,p,K,byrow=TRUE)),p,d)
  var.1=lambda*(1+disp*lambda)
  cl=(apply(f.z.y,1,which.max))
  time.sec= (proc.time() - ptm)[3]
  output<-list(y=y,K=K,cr=cr,cl=cl,likelihood=likelihood,AIC=AIC,BIC=BIC,a=a,
  lambda=lambda,f.z.y=f.z.y,w=w,time.sec=time.sec,variances=var.1,
  gname=gname)
  
  if (! quiet) {
    nomi.colonna= c("logL","BIC","AIC","EM-iter","Elapsed Time")
    df=data.frame(matrix(c(likelihood[length(likelihood)],round(BIC),round(AIC),h,time.sec),nrow=1))
    colnames(df)=nomi.colonna
    message("")
    message("Estimation details:")
    print(df)
    
  }
  
  invisible(output)
  
}


filter.em=function(y)
{
p=dim(y)[1]

matrice1=apply(y,1,sum)
matrice2=apply(y,1,mean)
cond=((matrice1>5) & (matrice2>0.5))
n.rm=length(which(!cond))
gname=(1:p)[cond]
out=list(n.rm=n.rm,y=y[cond,],gname=gname)
}




fun.a<-function(a,n,p,E.u.yz,E.logu.yz,f.z.y)
{
  preout1<-matrix(0,p,1)
  preout2<-matrix(0,p,1)
  
  preout1<-f.z.y*rowSums(E.logu.yz)
  preout2<-f.z.y*rowSums(E.u.yz)
  
  out<-n*(a*log(a)-lgamma(a))*sum(f.z.y)+(a-1)*sum(preout1)-a*sum(preout2)
  return(-out)
}


