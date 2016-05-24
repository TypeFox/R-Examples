book.stack.default<-
  function(x,B,k=2,alpha=0.05,bit=FALSE){
       
    options(warn=-1)
    check(test=3,x=x,B=B,k=k,alpha=alpha,bit=bit)
    
    res.tbl=book.stack.main(x=x,B=B,k=k,alpha=alpha,bit=bit)
    
    res.tbl$call = match.call()
    class(res.tbl) = c("book.stack","CryptRndTest")
    res.tbl
         
  }

