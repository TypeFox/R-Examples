latex.fdt_cat.multiple <- function(x,
                                   columns=1:6,
                                   where='!tbp',
                                   caption=NULL,
                                   label=NULL,
                                   size='', 
                                   algtable=c('\\flushleft', '\\centering', '\\flushright'),
                                   hline1='\\hline',
                                   header=NULL,
                                   hline2='\\hline',
                                   algclim='l',
                                   algfreq='r',
                                   hline3='\\hline')
{

 res <- list()
 for(i in 1:length(x)){
  y <- x[[i]][[1]]
  class(y) <- c('fdt_cat','data.frame') 
  res[[i]] <-  latex.fdt_cat(x=y,
                       columns=columns,
                       where=where,
                       caption=caption[i],
                       label=label[i],
                       size=size,
                       algtable=algtable,
                       hline1=hline1,
                       header=header,
                       hline2=hline2,
                       algclim=algclim,
                       algfreq=algfreq,
                       hline3=hline3) 

 }

 class(res) <- c('latex.fdt',
                 'list')

 return(res)
}
