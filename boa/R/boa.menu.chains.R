"boa.menu.chains" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nCHAINS MENU\n-----------"
   choices <- c("Back",
                "----------+",
                "Merge All |",
                "Delete    |",
                "Subset    |",
                "----------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = { if(boa.chain.collapse()) {
                    cat("+++ Successfully combined chains +++\n")
                 } else {
                    cat("--- No common parameters to combine ---\n")
                 }
               },
         "4" = { cat("\nDELETE CHAINS\n",
                     "=============\n",
                     "\nChains:\n",
                     "-------\n\n", sep = "")
                 lnames <- names(boa.chain("master"))
                 names(lnames) <- seq(lnames)
                 print(lnames)
                 value <- boa.getinput("\nSpecify chain index or vector of indices [none]\n")
                 if(length(value) > 0)
                    boa.chain.del(lnames = lnames[value])
               },
         "5" = { chain.args <- list()
                 info <- boa.chain.info(boa.chain("work"),
                                        boa.chain("work.support"))
                 cat("\nSUBSET CHAINS\n",
                     "=============\n",
                     "Specify the indices of the items to be included in the ",
                     "subset.  Alternatively,\n",
                     "items may be excluded by supplying negative indices.  ",
                     "Selections should be in\n",
                     "the form of a number or numeric vector.\n",
                     "\nChains:\n",
                     "-------\n\n", sep = "")
                 names(info$lnames) <- seq(info$lnames)
                 print(info$lnames)
                 value <- boa.getinput("\nSpecify chain indices [all]\n")
                 if(length(value) > 0) {
                    lnames <- info$lnames[value]
                    chain.args$lnames <- lnames
                 } else {
                    lnames <- info$lnames
                 }

                 cat("\nParameters:\n",
                     "-----------\n\n", sep = "")
                 info$pnames <- unique(unlist(info$pnames[lnames]))
                 names(info$pnames) <- seq(info$pnames)
                 print(info$pnames)
                 value <- boa.getinput("\nSpecify parameter indices [all]\n")
                 if(length(value) > 0)
                    chain.args$pnames <- info$pnames[value]

                 cat("\nIterations:\n",
                     "+++++++++++\n\n", sep = "")
                 print(info$iter.range[lnames, , drop = FALSE])
                 value <- boa.getinput("\nSpecify iterations [all]\n")
                 if(length(value) > 0)
                    chain.args$iter <- value

                 do.call("boa.chain.subset", args = chain.args)
               },
         "6" = NULL
      )
   }

   return(abs(idx))
}
