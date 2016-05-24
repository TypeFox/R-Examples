



.noinf <- function(df, clmns=names(df), buf=1){
  signs <- c(-1,1)
  min.max <- nv(c('min','max'), signs)

  for(clmn in clmns)
     for(sgn in signs){
       vct <- df[,clmn] 
       vct[!is.finite(vct) & sign(vct)==sgn]<- do.call(min.max[as.character(sgn)], list(vct[is.finite(vct)])) + sgn * buf
       df[,clmn] <- vct
     }
  df
}