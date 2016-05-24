##########################################################################
## this function automates the Scythe C++ call making book-keeping
## much easier
##
##    output.object:  name of posterior sample that will be placed
##                    in the parent environment (string)
##
##    cc.fun.name:    name of the C++ function to be called (string)
##
##    package:        name of package (string, "MCMCpack" by default)
##
##    developer:      option that determines whether the R call and
##                    C++ template are echoed or whether the Scythe
##                    call is made (logical)
##
##    help.file:      option that determines whether a template of a
##                    helpfile for the calling R function should be
##                    generated (logical)
##
##    cc.file:        the file used to store the C++ template
##                    (string, output to the screen if "")
##
##    R.file:         the file used to store the clean R template
##                    (string, output to the screen if "")
##                    this is just the function call indented nicely
##
##    ...:            list of objects passed to C++ 
##                    NOTE: this will only take integers (which have
##                    to be coerced), doubles, and matrices.  They
##                    should all be of the form "X = X," with the first
##                    part the C++ name and the second part the R name.
##                    Remember that C++ names cannot have periods in them.
##                    Matrices will be, for example, Xdata, Xrow, and
##                    Xcol.
##
##                    Any objects that are changed in the C++ code
##                    need to have C++ names like sample.nonconst.  All
##                    *.nonconst have the tag stripped and a pointer
##                    is used to change values.  This is used in all models
##                    for the posterior density sample (sample.nonconst),
##                    and other quantities of interest.
##
##       This also build a skeleton C++ template and clean R template
##       for MCMCpack if developer=TRUE.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Updated by ADM and KQ 1/25/2006 (to allow for multiple nonconsts)
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"auto.Scythe.call" <-
  function(output.object, cc.fun.name, package="MCMCpack",
           developer=FALSE, help.file=FALSE, cc.file="", R.file="", ...) {

   # pull stuff from function call
   objects <- list(...)
   K <- length(objects)
   c.names <- names(objects)
   nonconst.indic <- rep(FALSE, K)
   test <- grep("\\.nonconst$", c.names)
   if(length(test)==0)
      stop("something must be declared non-constant in auto.Scythe.call()\n")
   nonconst.indic[test] <- TRUE
   c.names <- sub("\\.nonconst$", "", c.names)
   if(length(unique(c.names)) != K)
      stop("non-unique nonconst names passed in auto.Scythe.call()\n")
      
   R.names <- strsplit(toString(match.call()), ",")[[1]]
   R.names <- R.names[(length(R.names)-K+1):length(R.names)]

   ## put default values in for burnin, mcmc, thin, and verbose
   ## if not explicitly supplied
   burnin.exist <- mcmc.exist <- thin.exist <- seed.exist <-
     verbose.exist <- FALSE
   for (k in 1:K){
     if(c.names[[k]]=="burnin" & is.integer(objects[[k]]))
       burnin.exist <- TRUE
     if(c.names[[k]]=="mcmc" & is.integer(objects[[k]]))
       mcmc.exist <- TRUE
     if(c.names[[k]]=="thin" & is.integer(objects[[k]]))
       thin.exist <- TRUE
     if(c.names[[k]]=="verbose" & is.null(dim(objects[[k]])))
       verbose.exist <- TRUE        
   }
   if (!burnin.exist){
     objects <- c(objects, burnin=as.integer(5000))
     R.names[length(objects)] <- "as.integer(5000)"
   }
   if (!mcmc.exist){
     objects <- c(objects, mcmc=as.integer(25000))
     R.names[length(objects)] <- "as.integer(25000)"     
   }
   if (!thin.exist){
     objects <- c(objects, thin=as.integer(1))
     R.names[length(objects)] <- "as.integer(1)"     
   }
   if (!verbose.exist){
     objects <- c(objects, verbose=as.integer(FALSE))
     R.names[length(objects)] <- "as.integer(FALSE)"     
   }

   ## check parameters
   check.mcmc.parameters(objects$burnin, objects$mcmc, objects$thin)
   
   ## write out a template R help file
   callfun <- strsplit(toString(sys.call(which=-1)),",")[[1]][1]
   if (help.file){
      prompt.call <- paste("prompt(", callfun, ", file =\"",
         paste(callfun, ".template.Rd\"", sep=""), ")")
      eval(parse(text=prompt.call))
   }
   
   ##
   ## pull together R call
   ##

   # strings for R call
   start <- paste(".C('", cc.fun.name, "',", sep="")
   end <- paste("PACKAGE='", package, "')", sep="") 
   middle <- NULL
   
   for(k in 1:K) {
     if(is.double(objects[[k]]) & is.null(dim(objects[[k]]))) {
       if (regexpr("as.double", R.names[[k]])==-1){
         holder <- paste(c.names[[k]], " = as.double(", R.names[[k]], "),",
                         sep="")
         middle <- c(middle, holder)
       }
       else {
         holder <- paste(c.names[[k]], " =", R.names[[k]], ",",
                         sep="")
         middle <- c(middle, holder)
       }
     }
     else if(is.integer(objects[[k]]) & is.null(dim(objects[[k]]))) {
       if (regexpr("as.integer", R.names[[k]])==-1){
         holder <- paste(c.names[[k]], " = as.integer(", R.names[[k]], "),",
                         sep="")
         middle <- c(middle, holder)
       }
       else {
         holder <- paste(c.names[[k]], " =", R.names[[k]], ",",
                         sep="")
         middle <- c(middle, holder)
       }
     }
     else if(is.matrix(objects[[k]])) {
       holder.data <- paste(c.names[[k]], "data", " = as.double(",
                            R.names[[k]], "),", sep="")
       holder.row <- paste(c.names[[k]], "row", " =", " nrow(",
                           R.names[[k]], "),", sep="")
       holder.col <- paste(c.names[[k]], "col", " =", " ncol(",
                           R.names[[k]], "),", sep="")
       middle <- c(middle, holder.data, holder.row, holder.col)
     }
     else {
       stop("Integers, doubles, or matrices only to auto.Scythe.call().")
     }
   }
   
   # clean up and return R call
   middle <- paste(middle, sep=" ", collapse=" ")
   call <- paste(start, middle, end, sep=" ")
   call <- gsub('\\( ', '\\(', call)  

   ##   
   ## pull together C++ call
   ##

   # strings for C++ call
   c.start <- paste("void ", cc.fun.name, "(", sep="")
   c.end <- ")"   
   c.middle <- NULL
   together.call <- NULL

   for(k in 1:K) {
      if(is.double(objects[[k]]) & is.null(dim(objects[[k]]))) {
         holder <- paste("double *", c.names[[k]], ",", sep="")
         if(!nonconst.indic[k]) holder <- paste("const ", holder, sep="")
         c.middle <- c(c.middle, holder)
      }
      else if(is.integer(objects[[k]]) & is.null(dim(objects[[k]]))) {
         holder <- paste("int *", c.names[[k]], ",", sep="")
         if(!nonconst.indic[k]) holder <- paste("const ", holder, sep="")         
         c.middle <- c(c.middle, holder)
      }
      else if(is.matrix(objects[[k]])) {
         holder.data <- paste("double *", c.names[[k]],
                               "data,", sep="")
         if(!nonconst.indic[k]) holder.data <-
            paste("const ", holder.data, sep="")
         scythe <- paste("     Matrix <double> ", c.names[[k]],
                         " = r2scythe(*",
                         c.names[[k]], "row, *", c.names[[k]], "col, ",
                         c.names[[k]], "data);\n", sep="")
                            
          together.call <- paste(together.call, scythe, sep="")
         
      holder.row <- paste("const int *", c.names[[k]], "row,", sep="")
      holder.col <- paste("const int *", c.names[[k]], "col,", sep="")
      c.middle <- c(c.middle, holder.data, holder.row, holder.col)
      }
      }

   # clean up and print C++ function call
   c.middle <- paste(c.middle, sep=" ", collapse=" ")
   c.call <- paste(c.start, c.middle, c.end, sep="")
   c.call <- gsub(',)', ')', c.call)

   # if developer dump Scythe code to file, R function to screen, and evaluate
   if(developer) {
      comment.block <- paste("// ", cc.file, " DESCRIPTION HERE\n//\n// The initial version of this file was generated by the\n// auto.Scythe.call() function in the MCMCpack R package\n// written by:\n//\n// Andrew D. Martin\n// Dept. of Political Science\n// Washington University in St. Louis\n// admartin@wustl.edu\n//\n// Kevin M. Quinn\n// Dept. of Government\n// Harvard University\n// kevin_quinn@harvard.edu\n// \n// This software is distributed under the terms of the GNU GENERAL\n// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE\n// file for more information.\n//\n// Copyright (C) ", substring(date(),21,24), " Andrew D. Martin and Kevin M. Quinn\n// \n// This file was initially generated on ", date(), "\n// REVISION HISTORY\n\n", sep="")
      
      includes.block <- '#include "matrix.h"\n#include "distributions.h"\n#include "stat.h"\n#include "la.h"\n#include "ide.h"\n#include "smath.h"\n#include "MCMCrng.h"\n#include "MCMCfcds.h"\n\n#include <R.h>           // needed to use Rprintf()\n#include <R_ext/Utils.h> // needed to allow user interrupts\n\nusing namespace SCYTHE;\nusing namespace std;\n\n'
      
      main.block <- 'extern "C" {\n\n   // BRIEF FUNCTION DESCRIPTION\n'
      
      function.call <- paste('   ', c.call, ' {\n', sep="")
      
      together.block <- "   \n     // pull together Matrix objects\n     // REMEMBER TO ACCESS PASSED ints AND doubles PROPERLY\n"

      constants.block <- "\n     // define constants\n     const int tot_iter = *burnin + *mcmc;  // total number of mcmc iterations\n     const int nstore = *mcmc / *thin;      // number of draws to store\n     const int NUMBER_OF_PARAMETERS = ????; // YOU NEED TO FILL THIS IN\n"

      storage.block <-   "\n     // storage matrix or matrices\n     Matrix<double> STORAGEMATRIX(nstore, NUMBER_OF_PARAMETERS);\n"      

      seed.block <- "\n     // initialize rng stream\n     rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);\n"

      startval.block <- "\n     // set starting values\n     PARAMETER_BLOCK1 = ????;\n     PARAMETER_BLOCK2 = ????;\n     ETC.;\n"
      
      sample.call <- paste('\n     ///// MCMC SAMPLING OCCURS IN THIS FOR LOOP\n     for(int iter = 0; iter < tot_iter; ++iter){\n\n       // sample the parameters\n       PARAMETER_BLOCK1 = ????;\n       PARAMETER_BLOCK2 = ????;\n       ETC;\n\n       // store draws in storage matrix (or matrices)\n       if(iter >= *burnin && (iter % *thin == 0)){\n        // PUT DRAWS IN STORAGEMATRIX HERE\n       }\n\n       // print output to stdout\n       if(*verbose == 1 && iter % 500 == 0){\n         Rprintf("\\n\\n', cc.fun.name, ' iteration %i of %i \\n", (iter+1), tot_iter);\n         // ADD ADDITIONAL OUTPUT HERE IF DESIRED\n       }\n\n       void R_CheckUserInterrupt(void); // allow user interrupts\n     } // end MCMC loop\n\n     delete stream; // clean up random number stream\n\n     // load draws into sample array\n', sep="")
      
      end.block <- paste("\n   } // end", cc.fun.name,  '\n} // end extern "C"\n')

      if (cc.file == ""){
        cat("\n\nThe C++ template file is:\n")
      }
      cat(comment.block, includes.block, main.block, function.call,
          together.block, together.call, constants.block, storage.block,
          seed.block, startval.block, sample.call,
          end.block, sep="", file=cc.file)
      if (cc.file != "") {
        cat("\nCreated file named '", cc.file, "'.\n", sep="")
        cat("Edit the file and move it to the appropriate directory.\n")
      }

      if (R.file == "") {
        cat("\n\nThe clean R template file is:\n")
      }
      dump(callfun, file=R.file)
      if (R.file != "") {
        cat("\nCreated file named '", R.file, "'.\n", sep="")
        cat("Edit the file and move it to the appropriate directory.\n")
        cat("Do not forget to edit the MCMCpack NAMESPACE file if\n")
        cat("installing new functions as part of MCMCpack.\n")
      }

      cat("\nThe call to .C is:\n")      
      draw.sample.call <- parse(text=paste(output.object, " <- ", call))
      print(draw.sample.call)
      cat("\nAUTOMATIC TEMPLATE FILE CREATION SUCCEEDED.\n")
      cat("All created templated files placed in ", getwd(), ".\n", sep="")
      invokeRestart("abort")
    }

   # if not developer evaluate call leaving output.object in
   # parent environment
   if(!developer) {
      draw.sample.call <- parse(text=paste(output.object, " <- ", call))
      eval(draw.sample.call, envir=parent.frame(1))
   }
}

