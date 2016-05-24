"boa.chain.del" <-
function(lnames, pnames)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   master <- boa.chain("master")
   master.support <- boa.chain("master.support")
   if(!missing(lnames)) {
      lnames <- intersect(lnames, names(master))
      for(i in lnames) {
         master[[i]] <- NULL
         master.support[[i]] <- NULL
      }
   }
   if(!missing(pnames)) {
      lnames <- names(master)
      for(i in lnames) {
         keep <- setdiff(boa.pnames(master[[i]]), pnames)
         master[[i]] <- master[[i]][, keep]
         master.support[[i]] <- master.support[[i]][, keep]
      }
   }
   if(boa.chain("work.sync")) {
      boa.chain(master = master, master.support = master.support,
                work = master, work.support = master.support)
   } else {
      boa.chain(master = master, master.support = master.support)
   }
   invisible()
}
