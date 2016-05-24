ptb <- function(x,y){ # x ist dichotom ; y ist kontinuierlich
  stopifnot(length(x)==length(y))
  index <- !is.na(x) & !is.na(y) # NAs rausnehmen (listwise)
  n=sum(index)
  x <- x[index]
  if((sum(range(x)==c(0,1)))==2 ){  
    # check ob x dichotom 0,1 und beides vertreten! 
    y <- y[index]
    M1 <- mean(y[x==1])
    M2 <- mean(y[x==0])
    p <- sum(x)/n
    q <- 1-p
    r <- ((M1-M2)*sqrt(p*q))/(sqrt(var(y)))
    return(list(rptb=r,n=n))
  }else{ return(list(rptb=NA,n=n)) }
}