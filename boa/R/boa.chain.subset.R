"boa.chain.subset" <-
function(lnames, pnames, iter)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   work.support <- boa.chain("work.support")
   if(!missing(lnames)) {
      for(i in setdiff(names(work), lnames)) {
         work[[i]] <- NULL
         work.support[[i]] <- NULL
      }
   }
   for(i in names(work)) {
      link <- work[[i]]
      link.support <- work.support[[i]]
      if(!missing(pnames)) {
         link <- boa.getparms(link, pnames)
         link.support <- boa.getparms(link.support, pnames)
      }
      if(!missing(iter))  link <- boa.getiter(link, iter)
      work[[i]] <- link
      work.support[[i]] <- link.support
   }
   subset <- length(work) > 0
   if(subset)  boa.chain(work = work, work.support = work.support,
                         work.sync = FALSE)

   return(subset)
}
