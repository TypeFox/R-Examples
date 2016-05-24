maxstat.setupSNP<-function(x, y, colSNPs=attr(x,"colSNPs"), ...)
{
 if(missing(y))
  stop("a case-control variable must be indicated in the second argument") 
 
 y<-deparse(substitute(y))
  if (!exists(y)) y<-x[,y] else y<-get(y)
 if(length(table(y))>2)
  stop("case-control variable must have only 2 levels")  

 ans<-sapply(x[,colSNPs, drop=FALSE], function(o) maxstat(y, o, ...))
 class(ans)<-"maxstat"
 ans
}
