###################################################
#### AUTHOR:     Arnost Komarek                ####
####             (2005)                        ####
####                                           ####
#### FILE:       bayessurvreg2.writeHeaders.R  ####
####                                           ####
#### FUNCTIONS:  bayessurvreg2.writeHeaders    ####
###################################################

### ======================================
### bayessurvreg2.writeHeaders
### ======================================
## Subfunction for bayessurvreg2.R
##  -> just to make it more readable
##
## Write headers to files where simulated values will be stored
##
bayessurvreg2.writeHeaders <- function(dir, doubly, prior.init, store, design, design2)
{
  bayesBisurvreg.writeHeaders(dir=dir, dim=1, nP=design$n, doubly=doubly,
                              prior.init=prior.init, store=store, design=design, design2=design2)

  FILES <- dir(dir)
  
  ## Files with sampled covariance matrices of random effects
  if (design$nrandom){
    sink(paste(dir, "/D.sim", sep = ""), append = FALSE)
    D <- diag(design$nrandom)
    rows <- row(D)[lower.tri(row(D), diag = TRUE)]
    cols <- col(D)[lower.tri(col(D), diag = TRUE)]            
    cat("det", paste("D.", rows, ".", cols, sep = ""), "\n", sep = "      ")
    sink() 
  }
  else{
    if ("D.sim" %in% FILES) file.remove(paste(dir, "/D.sim", sep = ""))
  }

  if (doubly){
    if (design2$nrandom){
      sink(paste(dir, "/D_2.sim", sep = ""), append = FALSE)
      D2 <- diag(design2$nrandom)
      rows2 <- row(D2)[lower.tri(row(D2), diag = TRUE)]
      cols2 <- col(D2)[lower.tri(col(D2), diag = TRUE)]            
      cat("det", paste("D.", rows2, ".", cols2, sep = ""), "\n", sep = "      ")
      sink() 
    }
    else{
      if ("D_2.sim" %in% FILES) file.remove(paste(dir, "/D_2.sim", sep = ""))
    }
  }
  else{
    if ("D_2.sim" %in% FILES) file.remove(paste(dir, "/D_2.sim", sep = ""))
  }    
  
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
}
