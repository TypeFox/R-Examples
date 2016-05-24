projsplx_k = function (y,K){
  rho_prin=0:K
  large_L=FALSE;
  N = length(y); bget = FALSE;
  if(K==N){
    return(rep(1/N,N))
  }
  s = sort(y,decreasing=TRUE); 
  y_sum=cumsum(s)
  ll=0
  
  while(ll<K+1){
    j_ind=(ll+1):N
    if(ll>0){
    y_sum_temp=y_sum[j_ind]-y_sum[ll]
    }else{
      y_sum_temp=y_sum
    }
    lambda = 1/(j_ind - ll) * (1 - ll/K - y_sum_temp);
    rho_check=s[j_ind] + lambda
    rho=tail(which(rho_check>0), n=1) 
    if(length(rho)>0){
      
        lambda = lambda[rho]
        if(length( which ( (abs(s[j_ind]+lambda) -1/K)>10^-8 )) > 0  ){
         # ind=which ( (abs(s[j_ind]+lambda) -1/K)>10^-8 )
        #  ttt=(abs(s[j_ind]+lambda) -1/K)
         # print(ttt[ind])
          ll=ll+1;
          next;
        }
        if(ll>0){
                  x_temp=s[1:ll]+lambda
                  if ( length(which( x_temp> (1/K) )) < ll){
                    ll=ll+1
                    next;
                  }
                  temp=(y+lambda)
                  temp[temp<0]=0
                  temp[temp>1/K]=1/K
                  x=temp;
                  return(x)
        } else   {
          rho_prin=0;
          temp=(y+lambda)
          temp[temp<0]=0
          x=temp;
          return(x)
        }

    }else{
      large_L=TRUE;
      return(c(rep(1/K,K),rep(0,N-K)))
    }
  
  }

  x=projsplx(y)
 return (x)
 
}