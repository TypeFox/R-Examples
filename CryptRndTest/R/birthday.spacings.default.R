birthday.spacings.default <-
  function(x,m=128,n=2^16,alpha=0.05,lambda,num.class=10){
    options(warn=-1)
    check(test=2,x=x,alpha=alpha,m=m,n=n,lambda=lambda,num.class=num.class)
     
    res.tbl=birthday.spacings.main(x,m=m,n=n,alpha=alpha,lambda,num.class=num.class)
    res.tbl$call = match.call()
    class(res.tbl) = c("birthday.spacings","CryptRndTest")
    res.tbl
    
  }