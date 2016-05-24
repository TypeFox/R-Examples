adaptive.chi.square.default <-
  function(x,B,S,alpha=0.05,prop=0.5,bit=FALSE){
    options(warn=-1)
    check(test=1,x=x,B=B,alpha=alpha,prop=prop,S=S,bit=bit)
    
    res.tbl=adaptive.chi.square.main(x=x,B=B,S=S,alpha=alpha,prop=prop,bit=bit)
    res.tbl$call = match.call()
    class(res.tbl) = c("adaptive.chi.square","CryptRndTest")                                          
    res.tbl
      
  }