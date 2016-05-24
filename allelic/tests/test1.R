library(allelic)

newallelicTable<- function (table) {
#  print(table);
  if ( ! is.matrix( table ) ) {
    if (length( table ) != 6 )
      stop("'table' must hold 6 integers, in row order");
    table <- matrix(table,nrow=2,ncol=3, byrow=TRUE)
  }
  if ( ! is.matrix( table ) || nrow(table) != 2 || ncol(table) != 3 )
    stop("'table' must be 2x3 matrix");
#  if (any(! is.integer(table) ) )
#    stop("'table' must hold integers");
  return( allelic.exact.test( table[1,1],table[1,2],table[1,3],table[2,1],table[2,2],table[2,3] ) );
}

systematicTestNewAllelic <- function(n1,n2,incr) {
#  print( c("n1",n1,"n2",n2,"incr",incr)  );
  table <- matrix(c(0,0,0,0,0,0),nrow=2,ncol=3, byrow=TRUE)
  pvalues <- c()
  fisher22 <- c()
  fisher23 <- c()
  table22 <- matrix(0,nrow=2,ncol=2, byrow=TRUE)
  for(a in seq(0,n1,by=incr)) {
    print(c("a=",a))
    table[1,1] <- a
    for(b in seq(0,n1-a,by=incr)) {
      table[1,2] <- b
      table22[1,1] <- a + a + b
      table[1,3] <- n1-a-b
      table22[1,2] <- table[1,3]*2 + b
      for(d in seq(0,n2,by=incr)) {
        table[2,1] <- d
       
        for(e in seq(0,n2-d,by=incr)) {
          table[2,2] <- e
          table[2,3] <- n2-d-e
          table22[2,1] <- d + d + e
          table22[2,2] <- table[2,3]*2 + e
          p <- newallelicTable(table)

          pvalues[length(pvalues)+1] <- p
        }
      }
    }
  }
  
  print( c("Got ",length(pvalues)," pvalues"))

  return( list(fisher22=fisher22,newallelic=pvalues,fisher23=fisher23))
}


l <- systematicTestNewAllelic(300,350,40)
sum <- sum(l$newallelic)
print(c("sum=",sum,"roundedSum=",round(sum)))

if (round(sum * 1000 ) != 92187)
  stop("bad value for the sum of pvalue, should have been rounded to 92.187")


