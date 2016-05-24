"boa.chain.eval" <- function(expr, pname)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   success <- TRUE
   
   master <- boa.chain("master")
   master.support <- boa.chain("master.support")
   for(i in names(master)) {
      x <- try(eval(expr, as.data.frame(master[[i]])), TRUE)
      if(inherits(x, "try-error")) {
         success <- FALSE
         cat(x)
         break
      }
      pnames <- boa.pnames(master[[i]])
      if(is.element(pname, pnames)) {
         master[[i]][, pname] <- x
         master.support[[i]][, pname] <- c(-Inf, Inf)
      } else {
         pnames <- c(pnames, pname)
         master[[i]] <- cbind(master[[i]], x)
         dimnames(master[[i]])[[2]] <- pnames
         master.support[[i]] <- cbind(master.support[[i]], c(-Inf, Inf))
         dimnames(master.support[[i]])[[2]] <- pnames
      }
   }
   if(boa.chain("work.sync")) {
      boa.chain(master = master, master.support = master.support,
                work = master, work.support = master.support)
   } else {
      boa.chain(master = master, master.support = master.support)
   }
   
   invisible(success)
}
