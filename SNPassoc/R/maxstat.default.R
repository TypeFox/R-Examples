maxstat.default<-function(x, y, ...)
 {
   # x: case/control
   # y: SNP

   name.caco<-deparse(substitute(x))
   name.snp<-deparse(substitute(y))
   ok <- complete.cases(x,y)
   x <- x[ok]
   y <- y[ok]

   if(length(unique(x))>2 ) stop(paste(name.caco, " has > 2 different values"))
   if(length(unique(y))>3 ) stop(paste(name.snp, " has > 3 different values"))
   xx<-table(x,y)
   maxstat.table(xx, ...)
 }
