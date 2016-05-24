kolmogorov<-function (x,fnc,...) 
{
  
  find.kol = function (d, n) 
  {
    sum.val=0
    sum.max = floor(n*(1-d))
    for(j in 0:sum.max){
      prod.1 = choose(n,j)
      prod.2 = (1-d-(j/n))^(n-j)
      prod.3 = (d+(j/n))^(j-1)
      sum.val = sum.val+(prod.1*prod.2*prod.3)
    }
    return(2*d*sum.val)
  }
  
  f.n = function(data,x){
    # F_{n}(x) returns # of X's in the sample(data) <= x / n
    return(length(data[data<=x])/length(data))  
  }
  x = sort(x)
  x.unique = sort(unique(x))
  
  if(length(x.unique)==length(x)){ #there are no ties
    
    D = 0
    n = length(x)
    f.0 = fnc(x,...)
    for(i in 1:n){
      D = max(D,abs((i/n)-f.0[i]),abs(((i-1)/n)-f.0[i])) 
    }
    
  }
  
  else{
    #there are ties
    
    fn = numeric()
    f.0 = fnc(x.unique,...)
    for(j in x.unique){
      fn = c(fn, f.n(x,j))
    }
    D = abs(fn[1]-f.0[1]) # i=1 
    for(i in 2:length(x.unique)){ 
      D = max(D,abs(fn[i]-f.0[i]),abs(fn[i-1]-f.0[i])) 
    }
  }
  
  p = find.kol(D,n)
  cat("D=", D, "\n", "p=", p,"\n")
  return(list(D=D,p=p))
  
}
