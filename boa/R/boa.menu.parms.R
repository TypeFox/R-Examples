"boa.menu.parms" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nPARAMETERS MENU\n---------------"
   choices <- c("Back",
                "-----------+",
                "Set Bounds |",
                "Delete     |",
                "New        |",
                "-----------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = { chain.args <- list()
                 info <- boa.chain.info(boa.chain("master"),
                                        boa.chain("master.support"))
                 cat("\nSET PARAMETER BOUNDS\n",
                     "====================\n",
                     "\nParameters:\n",
                     "-----------\n\n", sep = "")
                 info$pnames <- unique(unlist(info$pnames))
                 names(info$pnames) <- seq(info$pnames)
                 print(info$pnames)
                 value <- boa.getinput("\nSpecify parameter index or vector of indices [none]\n")
                 if(length(value) > 0) {
                    chain.args$pnames <- info$pnames[value]
                    value <- boa.getinput("\nSpecify lower and upper bounds as a vector [c(-Inf, Inf)]\n")
                    if(length(value) == 2) {
                       chain.args$limits <- value
                    } else {
                       chain.args$limits <- c(-Inf, Inf)
                    }
                    do.call("boa.chain.support", args = chain.args)
                 }
               },
         "4" = { info <- boa.chain.info(boa.chain("master"),
                                        boa.chain("master.support"))
                 cat("\nDELETE PARAMETERS\n",
                     "=================\n",
                     "\nParameters:\n",
                     "-----------\n\n", sep = "")
                 pnames <- unique(unlist(info$pnames))
                 names(pnames) <- seq(pnames)
                 print(pnames)
                 value <- boa.getinput("\nSpecify parameter index or vector of indices [none]\n")
                 boa.chain.del(pnames = pnames[value])
               },
         "5" = { info <- boa.chain.info(boa.chain("master"),
                                        boa.chain("master.support"))
                 cat("\nNEW PARAMETER\n",
                     "=============\n",
                     "\nCommon Parameters:\n",
                     "------------------\n\n", sep = "")
                 pnames <- info$pnames[[1]]
                 for(i in info$lnames[-1]) {
                    pnames <- intersect(info$pnames[[i]], pnames)
                 }
                 print(pnames)
                 cat("\nNew parameter name [none]\n")
                 value <- scan(what = "", n = 1, sep = "\n")
                 if(length(value) > 0) {
                    expr <- boa.getinput("\nDefine the new parameter as a function of the parameters listed above\n",
                                         evaluate = FALSE)
                    if(length(expr) > 0) boa.chain.eval(expr, value)
                 }
               },
         "6" = NULL
      )
   }

   return(abs(idx))
}
