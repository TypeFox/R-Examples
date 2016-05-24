compute.mc.t.statistic <-
function(x, y, mu=NULL, sigma){
  n = nrow(y)
  y.bar = apply (y, 2, mean)
  x1 = sum(x)
  
  sub.df = 0
  
  ## mu and beta
  if (is.null(mu)){ # estimating mu
    beta = (t(y)%*%x - y.bar*x1)/(sum(x^2) - x1^2/n)
    mu = y.bar - beta*x1/n
    
    se.adj = sqrt(sum((x-mean(x))^2)) # standard error adjustment when we estimate the intercept
  }else{ # no need to estimate mu
    beta = (t(y)%*%x - mu*x1)/sum(x^2)
    
    se.adj = sqrt(sum(x^2))
  }
  
  
  ## test statistic
  test.stat = beta/sigma*se.adj
  
  list (ts=abs(test.stat))
}
