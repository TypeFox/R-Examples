#################################################
###### 11-15-09: R package "LogitNet"
###### R functions for fitting LogitNet model 
#################################################

#################################################################################################
############## Derive the weight matrix for taking account of spatial correlations along genomic profiles in LogitNet model.


LogitNet.weight<- function(X.m, chr) 
{
 ############## some constant
  width=10  ### window size for loess smooth.
  
  lohdata=X.m
  lohch=chr

  N <- nrow(lohdata)  ### number of samples
  M <- ncol(lohdata)  ### number of genes

  chrom=unique(lohch)
  K <- length(chrom)

  # logistic regression for each locus
  logistic <- matrix(0, M, M)

  m.11=t(lohdata)%*%lohdata
  m.10=t(lohdata)%*%(1-lohdata)
  m.01=t(1-lohdata)%*%lohdata
  m.00=t(1-lohdata)%*%(1-lohdata)

  temp=m.11*m.00/m.10/m.01
  logistic=log(temp)
  logistic[is.na(logistic)]=0
  logistic[logistic==-Inf]=min(logistic[logistic!=-Inf])
  logistic[logistic==Inf]=max(logistic[logistic!=Inf])

##############
  weights <- matrix(0, nrow=M,ncol=M)

  for(chr in chrom)
  {
    K=sum(lohch==chr)
    cur.log=logistic[lohch==chr, lohch==chr]
    cur.result=matrix(0, K, K)

    for(i in 1:K)
       {
         cur=cur.log[i,]
         cur.sm=lowess(cur, f=min(1/3, width/K))$y 
         cut=median(abs(diff(cur.sm)))  
         if(i>1)
         { 
           left=1:(i-1)     
           if(sum(cur.sm[left]<cut)>0) 
            {   
              index=max(left[cur.sm[left]<cut])
              cur.sm[1:index]=0
            }
         } 
         if(i<K)
         {
            right=(i+1):K
            if(sum(cur.sm[right]<cut)>0) 
            {   
              index=min(right[cur.sm[right]<cut])
              cur.sm[index:K]=0
            }
          }
          cur.result[i,]=cur.sm
       }
      weights[lohch==chr, lohch==chr]=cur.result
   }
   w.s=(weights+t(weights))/2
   w.s=exp(w.s)   
   diag(w.s)=1
   return(w.s)
}

#################################################
############## Fit LogitNet for one lambda

LogitNet <- function(X.m, weight.m, lambda, beta.ini=NULL) 
{
 ################# some constant
  tol=0.1
  maxiter=500

 ################## 
  n <- nrow(X.m)
  p <- ncol(X.m)
  X.vector <- as.vector(t(X.m))
  weights<-weight.m*lambda
  weights.vector <- as.vector(weights)
  iter_output<- -1
  beta_output<-rep(0, p*p)

  if(is.null(beta.ini))
  {
    results<-.C("logistic_lasso_weighted_c",
              as.double(X.vector),
              as.integer(n),
              as.integer(p),
              as.double(weights.vector),
              as.double(tol),
              as.integer(maxiter),
              beta_output=as.double(beta_output),
              iter_output=as.integer(iter_output)
              )
  } else {
    beta0=as.vector(beta.ini)
    results<-.C("logistic_lasso_weighted_ini_c",
              as.double(X.vector),
              as.integer(n),
              as.integer(p),
              as.double(weights.vector),
              as.double(tol),
              as.integer(maxiter),
              as.double(beta0),
              beta_output=as.double(beta_output),
              iter_output=as.integer(iter_output)
              )
  }
  beta=matrix(results$beta_output, p, p, byrow=T)
  return(beta)
 }




####################################
############ Fit LogitNet with Cross-validation 

LogitNet.CV <- function(X.m, weight.m, lambda.v, fold =5) 
{
 ################# some constant
  tol=0.01
  maxiter=500
 ################## 

  n <- nrow(X.m)
  p <- ncol(X.m)
  l <- length(lambda.v)
  X.vector <- as.vector(t(X.m))
  likelihood.test <- matrix(0, fold, l)
  likelihood.train <- matrix(0, fold, l)
 
  beta_unbias <- array(0, dim=c(p,p,fold,l))
  beta_reg <- array(0, dim=c(p,p,fold,l))
  
  for(i in 1:l)
  {
    w<-weight.m*lambda.v[i]
    w.vector <- as.vector(w)

    beta1<-rep(0, p*p*fold)
    beta2<-rep(0, p*p*fold)
    like.test<-rep(0,fold)
    like.train<-rep(0,fold)

    results<-.C("crossvalid",
                as.double(X.vector),
                as.integer(n),
                as.integer(p),
                as.integer(fold),
                as.double(w.vector),
                as.double(tol),
                as.integer(maxiter),
                beta1=as.double(beta1),
                beta2=as.double(beta2),
                like.test=as.double(like.test),
                like.train=as.double(like.train)
                )
     
     likelihood.test[,i]=results$like.test
     likelihood.train[,i]=results$like.train
     beta_unbias[,,,i]=array(results$beta2, dim=c(p,p,fold))
     beta_reg[,,,i]=array(results$beta1, dim=c(p,p,fold))
  }
  list(beta_reg, beta_unbias, likelihood.test, likelihood.train)
}



