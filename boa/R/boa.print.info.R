"boa.print.info" <-
function(which = "work")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   switch(which,
      "work" = { chain <- boa.chain("work")
                 chain.support <- boa.chain("work.support")
                 cat("\n",
                     "WORKING CHAIN SUMMARY:\n",
                     "======================\n\n", sep = "") },
      "master" = { chain <- boa.chain("master")
                 chain.support <- boa.chain("master.support")
                 cat("\n",
                     "MASTER CHAIN SUMMARY:\n",
                     "=====================\n\n", sep = "") },
      chain <- NULL
   )

   chain.info <- boa.chain.info(chain, chain.support)
   if(is.list(chain.info)) {
      cat("Iterations:\n",
          "+++++++++++\n\n", sep = "")
      print(chain.info$iter.range)
      for(i in chain.info$lnames) {
         header <- paste("\nSupport: ", i, "\n", sep = "")
         cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
         print(chain.info$support[[i]])
      }
   } else {
      cat("Warning: chain contains no data\n")
   }
   invisible(chain.info)
}
