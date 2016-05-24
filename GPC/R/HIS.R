# HIS : Hyperbolic Index Set
# p=max degree order
# M=dimension

HIS <- function(p,M,q){
  A=c()
  A = rbind(A,matrix(c(rep(0,M-1),1),1,M))
  if(p>1){
    for (i in 2:p){
      A =rbind(A,c(rep(0,M-1),i));
      for (j in 2:min(i,M)){
        B = algo_H_Bis(i,j);
        C =cbind(matrix(0,length(B[,1]),M-j),B);
        A =rbind(A,C);
      }  
    }
  }
  A = t(apply(A,1,sort))
  Q_A=c()
  if(q<1){
    for( j in 1 : length(A[,1])){
      if(norm_q(A[j,],q) <= p){
        Q_A=rbind(Q_A,A[j,]);
      }
    }
  } else 
  { 
    Q_A=A
  }
    
  D=c()#algo_L(Q_A[1,]);
  for (i in 1:length(Q_A[,1])){    
    C=algo_L(Q_A[i,]);
    NumDiffValues=length(unique(Q_A[i,]))-1
    #if(length(C[,1]) != (choose(M,length(which(Q_A[i,]!=0)))*factorial(NumDiffValues))) print(i)
    D=rbind(D,C);
  }
  
  D = rbind(matrix(rep(0,M),1,M),D)
  return(D)

}

#################################

norm_q <- function(alpha,q){
alpha_q=(sum(alpha^q))^(1/q);
return(alpha_q)
}

#################################

algo_H<-function(p,j){
  a=rep(1,j)
  a[1] = p - j + 1;
  a[j+1] = -1;
  #A[1,]=a;
  A=a;
  
  while(a[2]<(a[1]-1)){
    a[1] = a[1] - 1;
    a[2] = a[2] + 1;
    k = 3;
    s = a[1] + a[2] - 1;
    while(a[k]>(a[1]-1)){
      s = s + a[k];
      k = k + 1;    
    }
    if( (k>j) ){
      A = rbind(A,a);
    }
    else{
      A = rbind(A,a);
      a2=a;
      x = a2[k] + 1 ;
      a2[k] = x;
      k = k - 1;
      while(k>1){
        a2[k] = x;
        s = s - x;
        k = k - 1;
      }
      a2[1] = s;
      #a2=a
      A = rbind(A,a2);
    }
  }
  
  if(length(dim(as.array(A)))==1){A=t(as.matrix(A));di=dim(A)} else{di=dim(A)}
  # alpha1=[zeros(di[1],1),A(:,1:di[2]-1)];
  if(di[1]==1){alpha1=t(as.matrix(A[,1:di[2]-1]))} else{alpha1=A[,1:di[2]-1]}
  alpha=alpha1
  for(i in 1 : di[1]){
    alpha[i,]=sort(alpha1[i,]);
  }
  # alpha=[];
  # for i = 1 : r
  #     B=algo_L( alpha1(i,:));
  #     alpha=[alpha;B];
  # end
  return(alpha)
}

algo_H_Bis<-function(p,j,a_In){
  if (nargs() == 2){
    a=rep(1,j)
  a[1] = p - j + 1;
  a=matrix(a,nrow=1)
  }
  else {
    a=a_In
    a=matrix(a,nrow=1)
  }
  #a[j+1] = -1;
  #A[1,]=a;
  A=a;
  A=matrix(A,nrow=1)
  while(a[2]<=(a[1]-1) ){
    #    if(a[2]==(a[1]-1) ) { if(length(a)>=3 & a[3]<=(a[2]-1)){a[1] = a[1] - 1;a[3] = a[3] + 1} else { return(A)}}
    if(a[2]==(a[1]-1)) {
      if(length(a)>=3)
      {
        l=0;
        while(l<=(length(a)-3) & a[3+l]==a[2]) {l= l+1}
        if(l < (length(a)-2) ) 
        {
          a[1] = a[1] - 1;
          a[3+l] = a[3+l] + 1
        }
        else {
          return(A)
        }
      }
      else { 
        return(A)
      }
    }
    else {
      a[1] = a[1] - 1;
      a[2] = a[2] + 1;
    }
    #resulTemp = algo_H_Bis(p-a[1],j-1,a[-1])
    if (j>2){
      resulTemp = algo_H_Bis(p-a[1],j-1)
      if(any(resulTemp[,1]<=a[1]) == TRUE)
      {
        IndexTmp=which(resulTemp[,1]<=a[1])
        resulTemp=matrix(resulTemp[ IndexTmp,],nrow=length(IndexTmp))
        a=cbind(matrix(rep(a[1],dim(resulTemp)[1]),ncol=1),resulTemp)
      }
    }
    A = rbind(A,a)
    #a=a[nrow(resulTemp),]    
    #if (nrow(a)>1) {for(i in 1:nrow(a)) {Atemp = algo_H_Bis(p,j,a[i,]) ; A = rbind(A,Atemp)}}
    if (j>2) a=matrix(a[nrow(resulTemp),],nrow=1)  

  }
  return(A)
}

#################################
algo_L <- function(alpha){
  #di=dim(alpha); #di1=h; di2=c
  m=length(alpha);
  k=m-1;
  all_alpha=alpha;
  
  while(k>0){
    k
    while(alpha[k]>=alpha[k+1]){
      k=k-1;
      if(k==0){
        return(all_alpha)
      }
    }
    l=m;
    while(alpha[k]>=alpha[l]){
      l=l-1;
    }
    an1=alpha[k];
    an2=alpha[l];
    alpha[k]=an2;
    alpha[l]=an1;
    r=k+1;
    l=m;
    while(r<l){
      an1=alpha[r];
      an2=alpha[l];
      alpha[r]=an2;
      alpha[l]=an1;
      r=r+1;
      l=l-1;
    }
    all_alpha=rbind(all_alpha,alpha);
    k=m-1;
  }
  return(all_alpha)
}