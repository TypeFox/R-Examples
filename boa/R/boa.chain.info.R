"boa.chain.info" <-
function(chain, chain.support)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   result <- NULL
   if(length(chain) > 0) {
      lnames <- names(chain)
      pnames <- list()
      iter <- list()
      iter.range <- NULL
      support <- list()
      for(i in lnames) {
         pnames[[i]] <- boa.pnames(chain[[i]])
         iter[[i]] <- boa.iter(chain[[i]])
         iter.range <- rbind(iter.range, c(range(iter[[i]]), length(iter[[i]])))
         support[[i]] <- chain.support[[i]]
      }
      dimnames(iter.range) <- list(lnames, c("Min", "Max", "Sample"))
      result <- list(lnames = lnames, pnames = pnames, iter = iter,
                     iter.range = iter.range, support = support)
   }

   return(result)
}
