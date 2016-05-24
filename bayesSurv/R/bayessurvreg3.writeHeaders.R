####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2005)                         ####
####                                            ####
#### FILE:       bayessurvreg3.writeHeaders.R   ####
####                                            ####
#### FUNCTIONS:  bayessurvreg3.writeHeaders     ####
####################################################

### ======================================
### bayessurvreg3.writeHeaders
### ======================================
## Subfunction for bayessurvreg3.R
##  -> just to make it more readable
##
## Write headers to files where simulated values will be stored
##
bayessurvreg3.writeHeaders <- function(dir, doubly, prior.init, priorb.di, priorb2.di, store, design, design2, version, mclass)
{
  bayesBisurvreg.writeHeaders(dir=dir, dim=1, nP=design$n, doubly=doubly,
                              prior.init=prior.init, store=store, design=design, design2=design2)

  FILES <- dir(dir)
  
  ## Files with sampled values of random effects
  if (store$b){
    sink(paste(dir, "/b.sim", sep = ""), append = FALSE)
    cat(paste(rep(design$names.random, design$ncluster), ".",
              rep(unique(design$cluster), rep(design$nrandom, design$ncluster)), sep = ""), "\n", sep = "    ")
    sink()
  }
  else{
    if ("b.sim" %in% FILES) file.remove(paste(dir, "/b.sim", sep = ""))
  }

  if (doubly){
    if (store$b2){
      sink(paste(dir, "/b_2.sim", sep = ""), append = FALSE)
      cat(paste(rep(design2$names.random, design2$ncluster), ".",
                rep(unique(design2$cluster), rep(design2$nrandom, design2$ncluster)), sep = ""), "\n", sep = "    ")
      sink()
    }
    else{
      if ("b_2.sim" %in% FILES) file.remove(paste(dir, "/b_2.sim", sep = ""))
    }    
  }
  else{
    if ("b_2.sim" %in% FILES) file.remove(paste(dir, "/b_2.sim", sep = ""))
  }

  ## Files with sampled G-splines - distribution of the random intercept
  if (version == 32){
    clean.Gspline(dir=dir, label="_b", care.of.y=FALSE)
    clean.Gspline(dir=dir, label="_b2", care.of.y=FALSE)

    sink(paste(dir, "/D.sim", sep = ""), append = FALSE)
    D <- diag(2)
    rows <- row(D)[lower.tri(row(D), diag = TRUE)]
    cols <- col(D)[lower.tri(col(D), diag = TRUE)]            
    cat("det", paste("D.", rows, ".", cols, sep = ""), "\n", sep = "      ")
    sink()         
  }
  else{
    if ("D.sim" %in% FILES) file.remove(paste(dir, "/D.sim", sep = ""))
    
    if (design$nrandom){
      write.headers.Gspline(dir=dir, dim=1, nP=design$ncluster, label="_b", gparmi=priorb.di$GsplI,
                            store.a=store$a.b, store.y=FALSE, store.r=store$r.b, care.of.y=FALSE)
    }
    else{
      clean.Gspline(dir=dir, label="_b", care.of.y=FALSE)
    }    

    if (doubly){
      if (design$nrandom){
        write.headers.Gspline(dir=dir, dim=1, nP=design2$ncluster, label="_b2", gparmi=priorb2.di$GsplI,
                              store.a=store$a.b2, store.y=FALSE, store.r=store$r.b2, care.of.y=FALSE)
      }
      else{
        clean.Gspline(dir=dir, label="_b2", care.of.y=FALSE)
      }     
    }
    else{
      clean.Gspline(dir=dir, label="_b2", care.of.y=FALSE)
    }
  }  

  ## Files with sampled values of correlation coefficient between the random intercepts and its Fisher Z transform
  if (version == 31){
    sink(paste(dir, "/rho_b.sim", sep = ""), append = FALSE)
    cat("rho    Z\n")
    sink()
  }
  else{
    if ("rho_b.sim" %in% FILES) file.remove(paste(dir, "/rho_b.sim", sep = ""))
  }  

  ## Files related to misclassification model
  if (mclass$nModel > 0){
    sink(paste(dir, "/sens_spec.sim", sep = ""), append = FALSE)    
    switch (mclass$Model,
      "Examiner" = {
        cat(paste(paste("alpha", mclass$labelExaminer, sep = ""), collapse = "  "))
        cat("  ")
        cat(paste(paste("eta",   mclass$labelExaminer, sep = ""), collapse = "  "))
        cat("\n")
      },
      "Factor:Examiner" = {
        cat(paste(paste("alpha", rep(mclass$labelExaminer, each = mclass$nFactor), ".", rep(mclass$labelFactor, mclass$nExaminer), sep = ""), collapse = "  "))
        cat("  ")
        cat(paste(paste("eta", rep(mclass$labelExaminer, each = mclass$nFactor), ".", rep(mclass$labelFactor, mclass$nExaminer), sep = ""), collapse = "  "))        
        cat("\n")
      },
      {
        sink()
        stop("some part of the code not implemented for this option")
      }  
    )
    sink()

    sink(paste(dir, "/devianceGJK.sim", sep = ""), append = FALSE)
    cat("D\n")
    sink()
  }  
}
