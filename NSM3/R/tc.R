tc<-function (x, tau, alternative="two.sided") 
{
  # A function for section 11.4: A Test of Exponentiality versus 
  # a trend change in mean residual life (Guess-Hollander-Proschan)
  
  f.n = function(data,x){
    # F_{n}(x) returns # of X's in the sample(data) <= x / n
    return(length(data[data<=x])/length(data))  
  }
  
  b.1 = function(u, data, tau){
    f.n.tau = f.n(data=data, x=tau)
    a = (2/3) - f.n.tau + ((f.n.tau^2)/2)
    b = -1 + f.n.tau - ((f.n.tau^2)/2)
    ans = (a*u) + (b*(u^2)) + ((u^4)/3)
    return(ans)
  }
  
  b.2 = function(u, data, tau){
    f.n.tau = f.n(data=data, x=tau)  
    a = (-1/6) + (f.n.tau/2) - ((f.n.tau^2)/2) + ((f.n.tau^3)/3) 
    b = (1/2) - f.n.tau + ((f.n.tau^2)/2)
    ans = (a*u) + (b*(u^2)) - ((u^4)/3)
    return(ans)
  }
  
  find.T = function(x, tau){
    # step 1, figure out if there are ties
    x.tilda = unique(x)
    x.sort = sort(x)
    n = length(x)
    if(length(x.tilda)==length(x)){ # There are no ties
      i.star = length(x.sort[x.sort<=tau])
      
      sum.1 = b.1(1, x.sort, tau)*x.sort[1] # i = 1
      for(i in 2:i.star){
        sum.1 = sum.1 + (b.1((n-i+1)/n, x.sort,tau)*(x.sort[i]-x.sort[i-1]))
      }
      
      sum.2 = 0
      for(i in (i.star+2):n){
        sum.2 = sum.2 + (b.2(((n-i+1)/n),x.sort,tau)*(x.sort[i]-x.sort[i-1]))
      }
      
      T1 = sum.1 + (b.1((n-i.star)/n, x.sort, tau)*(tau-x.sort[i.star])) + 
        (b.2((n-i.star)/n, x.sort, tau)*(x.sort[i.star+1]-tau)) + sum.2
      
      return(T1)
    }
    
    else{ # There are ties  
      x.tilda.sort = sort(x.tilda)
      i.star = length(x.tilda.sort[x.tilda.sort<=tau])
      k = length(x.tilda)
      # calculate all of the s_{i}'s once
      s.i = numeric()
      for(i in 1:k){
        idx = which(x.sort==x.tilda.sort[i])
        s.i = c(s.i, length(x.sort[idx]))
      }
      s.i = n-cumsum(s.i)
      sum.1 = b.1((s.i[1]-1)/n,x.tilda.sort, tau)*x.tilda.sort[1] # i=1
      for(i in 2:i.star){
        sum.1 = sum.1 + (b.1((s.i[i]-1)/n,x.tilda.sort, tau)*
                           (x.tilda.sort[i]-x.tilda.sort[i-1]))
      }
      sum.2 = 0
      for(i in (i.star +2):k){
        sum.2 = sum.2 + (b.2((s.i[i]-1)/n, x.tilda.sort, tau)*
                           (x.tilda.sort[i]-x.tilda.sort[i-1]))
      }
      T1 = sum.1 + (b.1(s.i[i.star]/n,x.tilda.sort, tau)*(tau-x.tilda.sort[i.star]))+
        (b.2(s.i[i.star]/n, x.tilda.sort,tau)*(x.tilda.sort[i.star+1]-tau)) + sum.2
      return(T1)
    }
  }
  
  
  
  
  # The next three lines are modified from fisher.test to get the correct alternative
  # hypotheses.  Citation needed?
  alternative <- char.expand(alternative, c("two.sided", "idmrl", "dimrl"))
  if (length(alternative) > 1L || is.na(alternative)) 
    stop("alternative must be \"two.sided\", \"idmrl\" or \"dimrl\"")
  
  T1 = find.T(x, tau)
  n = length(x)
  
  # large sample approximation
  
  f.n.tau = f.n(x, tau)
  sigma.sq = (mean(x)^2)* ( (-(f.n.tau^5)/15) + ((f.n.tau^4)/6) - ((f.n.tau^3)/6) +
                              ((f.n.tau^2)/10) - (f.n.tau/30) + (1/210) )
  sigma = sqrt(sigma.sq)
  T1.star=(sqrt(n)*T1/sigma)
  #dimrl
  if(alternative=="dimrl"){
    p=pnorm(T1.star)
  }
  
  #idmrl
  if(alternative=="idmrl"){
    p=pnorm(T1.star, lower.tail=F)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*pnorm(abs(T1.star), lower.tail=F)
  }
  
  cat("T1=", T1, "\n", "p=", p,"\n" , "sigma hat=", sigma,"\n" , "T1*=",T1.star ,"\n")
  return(list(T1=T1,prob=p, sigma.hat=sigma, T1.star=T1.star))
  
}
