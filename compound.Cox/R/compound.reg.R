compound.reg <-
function(t.vec,d.vec,X.mat,K=5,delta_a=0.025,a_0=0,var=FALSE,plot=TRUE,randomize=TRUE,var.detail=FALSE){

  n=length(t.vec)
  p=ncol(X.mat)
  
  if(randomize==TRUE){
    rand=sample(1:n)  ### randomize the subject IDs for cross-validation ###
    t.vec=t.vec[rand]
    d.vec=d.vec[rand]
    X.mat=X.mat[rand,]
  } 

  t.ot=t.vec[order(t.vec)]
  atr_t=(matrix(t.vec,n,n,byrow=TRUE)>=matrix(t.vec,n,n))
  
  l.func=function(b){
    X.matb=X.mat%*%b
    l1=sum( (X.matb)[d.vec==1] )
    S0=sum( (log(atr_t%*%exp(X.matb)))[d.vec==1] )
    l1=l1-S0
    as.numeric( l1 )
  }
  
  ###### Define Cross validated likelihood #######
  CV_a.func=function(a){
    CV_a=0
    
    for(k in 1:K){
      temp=(1+(k-1)*n/K):(k*n/K) ### index for the kth fold ###
      t.cv=t.vec[-temp];d.cv=d.vec[-temp];X.cv=X.mat[-temp,]
      m=length(t.cv)
      atr=(matrix(t.cv,m,m,byrow=TRUE)>=matrix(t.cv,m,m))
      
      la.func_cv=function(b){
        X.cvb=X.cv%*%b   
        l1=l0=sum( (X.cvb)[d.cv==1] )
        S0=sum( (log(atr%*%exp(X.cvb)))[d.cv==1] ) 
        l1=l1-S0
        bmat=matrix(b,m,p,byrow=TRUE)
        X.cvb0=X.cv*bmat
        S0_vec=colSums((log(atr%*%exp(X.cvb0)))*matrix(d.cv,m,p))
        l0=l0-sum( S0_vec ) 
        as.numeric(-( a*l1+(1-a)*l0 ))
      }
      
      beta_a_cv=nlm(la.func_cv,p=rep(0,p))$estimate
      
      X.cvba=X.cv%*%beta_a_cv
      l1=sum( (X.cvba)[d.cv==1] )
      S0=sum(  (log(atr%*%exp(X.cvba)))[d.cv==1]  )
      l_cv=l1-S0
      
      CV_a=CV_a+l.func(beta_a_cv)-l_cv
    }
    
    CV_a
  }
  
  ### Grid search on CV ####
  a_old=a_0
  CV_old=CV_a.func(a_old)
  CV_stock=CV_old
  repeat{
    a_new=a_old+delta_a
    CV_new=CV_a.func(a_new)
    CV_stock=c(CV_stock,CV_new)
    if(CV_new<=CV_old){break}
    CV_old=CV_new;a_old=a_new
    if(a_old>=1){break}
  }
  a_hat=min(a_old,1)

  if(plot==TRUE){
    plot(seq(a_0,a_new,length=length(CV_stock)),CV_stock,xlab="a",ylab="CV(a)",type="b")
    points(a_hat,CV_old,col="red",pch=17)
    abline(v=a_hat,lty="dotted")
  }

  ###### Shrinkage analysis ########
  la.func=function(b){
    X.matb=X.mat%*%b   
    l1=l0=sum( (X.matb)[d.vec==1] )
    S0=sum( (log(atr_t%*%exp(X.matb)))[d.vec==1] ) 
    l1=l1-S0
    bmat=matrix(b,n,p,byrow=TRUE)
    X.matb0=X.mat*bmat
    S0_vec=colSums((log(atr_t%*%exp(X.matb0)))*matrix(d.vec,n,p))
    l0=l0-sum( S0_vec ) 
    as.numeric(-( a_hat*l1+(1-a_hat)*l0 ))
  }
  
  res_a=nlm(la.func,p=rep(0,p),hessian=var)
  beta_a=res_a$estimate
  
  names(beta_a)=colnames(X.mat)
 
  ####### Variance estimation ########
  if(var==TRUE){
    V_hat=res_a$hessian/n
    b=beta_a
    h_dot=rep(0,p)
    for(i in 1:n){
      if(d.vec[i]==1){
        temp=t.vec>=t.vec[i]
        X_at=matrix(X.mat[temp,],nrow=sum(temp))
        S0=sum(  exp(X_at%*%b)  )
        S1=t(X_at)%*%exp(X_at%*%b)
        b_mat=matrix(b,sum(temp),p,byrow=TRUE)
        S0_vec=colSums(exp(X_at*b_mat))
        S1_vec=colSums(X_at*exp(X_at*b_mat))
        h_dot=h_dot-S1/S0+S1_vec/S0_vec
      } 
    }
    h_dot=h_dot/n
  
    dd_CV=as.numeric( hessian(CV_a.func,a_hat) )
    E_z=-dd_CV/n
    v_CV=t(h_dot)%*%solve(V_hat)%*%h_dot/(E_z^2)
    A_hat=solve(V_hat)%*%(  h_dot%*%t(h_dot)  )/E_z+diag(p)
    Sigma_hat=A_hat%*%solve(V_hat)%*%t(A_hat)
  
    SD_beta=sqrt(  diag(Sigma_hat)/n  )
    LCI=beta_a-1.96*SD_beta
    UCI=beta_a+1.96*SD_beta
    
    if(var.detail==TRUE){
      list(a_hat=a_hat,beta_hat=beta_a,SD=SD_beta,Lower95CI=LCI,Upper95CI=UCI,
           Sigma_hat=Sigma_hat,V_hat=V_hat,Hessian_CV=dd_CV,h_dot=as.vector(h_dot))
    }
    else{list(a_hat=a_hat,beta_hat=beta_a,SD=SD_beta,Lower95CI=LCI,Upper95CI=UCI)}
  }
  else{list(a_hat=a_hat,beta_hat=beta_a)}
  
}
