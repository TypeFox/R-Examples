"boa.chain.support" <-
function(lnames, pnames, limits)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   changed.master <- FALSE
   changed.work <- FALSE
   master <- boa.chain("master")
   master.support <- boa.chain("master.support")
   if(missing(lnames)) {
      lnames.master <- names(master.support)
   } else {
      lnames.master <- intersect(names(master.support), lnames)
   }
   for(i in lnames.master) {
      ipnames <- intersect(boa.pnames(master.support[[i]]), pnames)
      for(j in ipnames) {
         prange <- range(master[[i]][, j])
         if((limits[1] <= prange[1]) && (prange[2] <= limits[2])) {
            changed.master <- TRUE
            master.support[[i]][, j] <- limits
         } else {
            cat("Warning: support for '", j, "' in '", i, "' must include (",
                prange[1], ", ", prange[2], ")\n", sep ="")
         }
      }
   }
   work.sync <- boa.chain("work.sync")
   if(changed.master && work.sync) {
      changed.work <- TRUE
      work.support <- master.support
   } else if(!work.sync) {
      work <- boa.chain("work")
      work.support <- boa.chain("work.support")
      lnames.work <- names(work.support)
      if(!missing(lnames)) {
         idx <- pmatch(paste(lnames, "::", sep = ""),
                     paste(lnames.work, "::", sep = ""), nomatch = 0)
         lnames.work <- lnames.work[sort(idx)]
      }
      for(i in lnames.work) {
         ipnames <- intersect(boa.pnames(work.support[[i]]), pnames)
         for(j in ipnames) {
            prange <- range(work[[i]][, j])
            if((limits[1] <= prange[1]) && (prange[2] <= limits[2])) {
               changed.work <- TRUE
               work.support[[i]][, j] <- limits
            } else {
               cat("Warning: support for '", j, "' in 'work$", i,
                   "' must include (", prange[1], ", ", prange[2], ")\n",
                   sep ="")
            }
         }
      }
   }
   if(changed.master)  boa.chain(master.support = master.support)
   if(changed.work)  boa.chain(work.support = work.support)
   changed <- c(changed.master, changed.work)
   names(changed) <- c("master", "work")

   return(changed)
}
