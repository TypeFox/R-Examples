####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2005)                         ####
####                                            ####
#### FILE:       bayesBisurvreg.writeHeaders.R  ####
####                                            ####
#### FUNCTIONS:  bayesBisurvreg.writeHeaders    ####
####################################################

### ======================================
### bayesBisurvreg.writeHeaders
### ======================================
## Subfunction for bayesBisurvreg.R
##  -> just to make it more readable
##
## Write headers to files where simulated values will be stored
##
bayesBisurvreg.writeHeaders <- function(dir, dim, nP, doubly, prior.init, store, design, design2)
{
  FILES <- dir(dir)
  
  ## Common files
  sink(paste(dir, "/iteration.sim", sep = ""), append = FALSE)
  cat("iteration", "\n"); sink()

  ## Files for the G-spline related quantities   
  write.headers.Gspline(dir=dir, dim=dim, nP=nP, label="", gparmi=prior.init$Gparmi, store$a, store$y, store$r, care.of.y=TRUE)
  if (doubly) write.headers.Gspline(dir=dir, dim=dim, nP=nP, label="_2", gparmi=prior.init$Gparmi2,
                                    store$a2, store$y2, store$r2, care.of.y=TRUE)
  else        clean.Gspline(dir=dir, label="_2", care.of.y=TRUE)

  ## Files for regression parameters beta
  if (design$nX){ sink(paste(dir, "/beta.sim", sep = ""), append = FALSE)
                  cat(colnames(design$X), "\n", sep = "      "); sink() }
  else{
    if ("beta.sim" %in% FILES) file.remove(paste(dir, "/beta.sim", sep = ""))
  }

  if (doubly){
    if (design2$nX){ sink(paste(dir, "/beta_2.sim", sep = ""), append = FALSE)
                     cat(colnames(design2$X), "\n", sep = "      "); sink() }
    else{
      if ("beta_2.sim" %in% FILES) file.remove(paste(dir, "/beta_2.sim", sep = ""))
    }
  }
  else{
    if ("beta_2.sim" %in% FILES) file.remove(paste(dir, "/beta_2.sim", sep = ""))
  }  
}  
