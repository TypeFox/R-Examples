test_info <- function(ermobject,theta=seq(-5,5,0.01))
##Calculates info of a scale of items
#
#@input: ermobject ... Object of class eRm
#        theta ... supporting or sampling points on latent trait
#@output: a list where each element corresponds to an item and contains
#         $c.info...matrix of category information in columns for theta (rows)  
#         $i.info...vector of item information at values of theta
#@author: Thomas Rusch
#@date:12.6.2011
#  
  {
   infos <- item_info(ermobject,theta)
   tmp <- lapply(infos, function(x) x$i.info)
   tmp <- matrix(unlist(tmp),ncol=dim(ermobject$X)[2])
   tinfo <- rowSums(tmp)
   return(tinfo)
 }
