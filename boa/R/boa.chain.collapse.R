"boa.chain.collapse" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   lnames <- names(work)
   pnames <- boa.pnames(work[[1]])
   for(i in lnames[-1]) {
      pnames <- intersect(boa.pnames(work[[i]]), pnames)
   }
   link <- NULL
   for(i in lnames) {
      link <- rbind(link, boa.getparms(work[[i]], pnames))
   }
   collapse <- is.matrix(link)
   if(collapse) {
      lname <- paste(lnames, collapse = "::")
      work <- list(boa.sortiter(link))
      names(work) <- lname
      work.support <- list(boa.getparms(boa.chain("work.support")[[1]],
                                        boa.pnames(work[[1]])))
      names(work.support) <- lname
      boa.chain(work = work, work.support = work.support, work.sync = FALSE)
   }

   return(collapse)
}
