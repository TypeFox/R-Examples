# Produces a summary table for ps object 
summary.ps <- function(object,...){
      summary.tab <- NULL
      typ <- NULL   
   
      n.tp <- length(object$desc)
      for(i.tp in 1:n.tp){
         desc.temp <- object$desc[[i.tp]]
         iter      <- desc.temp$n.trees
         tp        <- names(object$desc)[i.tp]

		summary.tab <- rbind(summary.tab,
            with(desc.temp, c(n.treat,n.ctrl,ess.treat,
                                       ess.ctrl,
                                       max.es,
                                       mean.es,
                                       max.ks,
                                       max.ks.p,
                                       mean.ks,
                                       iter)))

                                       
                                       typ <- c(typ, tp)
      }
      

summary.tab <- matrix(summary.tab, nrow = n.tp)
      rownames(summary.tab) <- typ
      colnames(summary.tab) <- c("n.treat", "n.ctrl", "ess.treat", "ess.ctrl", "max.es", "mean.es", "max.ks", "max.ks.p","mean.ks","iter")
      
       
      class(summary.tab) <- "summary.ps"
      return(summary.tab)
}