####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2004)                         ####
####                                            ####
#### FILE:       bayessurvreg1.writeHeaders.R   ####
####                                            ####
#### FUNCTIONS:  bayessurvreg1.writeHeaders     ####
####################################################

### ======================================
### bayessurvreg1.writeHeaders
### ======================================
## Subfunction for bayessurvreg1.R
##  -> just to make it more readable
##
## Write headers to files where simulated values will be stored
##
bayessurvreg1.writeHeaders <- function(dir, prior, store, nX, X, names.random, ncluster, nrandom,
                                       rnamesX, unique.cluster, nBetaBlocks, nbBlocks)
{
   FILES <- dir(dir)
  
   sink(paste(dir, "/iteration.sim", sep = ""), append = FALSE)
   cat("iteration", "\n"); sink()

   sink(paste(dir, "/loglik.sim", sep = ""), append = FALSE)
   cat("loglik", "randomloglik", "\n", sep = "  "); sink()

   sink(paste(dir, "/mweight.sim", sep = ""), append = FALSE)
   cat(paste("w", 1:prior$kmax, sep = ""), "\n", sep = "      "); sink()

   sink(paste(dir, "/mmean.sim", sep = ""), append = FALSE)
   cat(paste("mu", 1:prior$kmax, sep = ""), "\n", sep = "      "); sink()

   sink(paste(dir, "/mvariance.sim", sep = ""), append = FALSE)
   cat(paste("sigma2", 1:prior$kmax, sep = ""), "\n", sep = "      "); sink()   
   
   sink(paste(dir, "/mixmoment.sim", sep = ""), append = FALSE)
   cat("       k", "              Intercept", "           Scale", "\n", sep = "  "); sink()
   
   if (nX){ sink(paste(dir, "/beta.sim", sep = ""), append = FALSE)
            cat(colnames(X), "\n", sep = "      "); sink() }
   else
     if ("beta.sim" %in% FILES) file.remove(paste(dir, "/beta.sim", sep = ""))
   

   if (store$b){ sink(paste(dir, "/b.sim", sep = ""), append = FALSE)
                 cat(paste(rep(names.random, ncluster), ".", rep(unique.cluster, rep(nrandom, ncluster)), sep = ""),
                     "\n", sep = "    "); sink() }
   else
     if ("b.sim" %in% FILES) file.remove(paste(dir, "/b.sim", sep = ""))

   if (nrandom){ sink(paste(dir, "/D.sim", sep = ""), append = FALSE)
                 D <- diag(nrandom)
                 rows <- row(D)[lower.tri(row(D), diag = TRUE)]
                 cols <- col(D)[lower.tri(col(D), diag = TRUE)]            
                 cat("det", paste("D.", rows, ".", cols, sep = ""), "\n", sep = "      "); sink() }
   else
     if ("D.sim" %in% FILES) file.remove(paste(dir, "/D.sim", sep = ""))

   if (store$y){sink(paste(dir, "/Y.sim", sep = ""), append = FALSE)
                cat(paste("Y", rnamesX, sep = ""), "\n", sep = "      "); sink() }
   else
     if ("Y.sim" %in% FILES) file.remove(paste(dir, "/Y.sim", sep = ""))   

   if (store$r){sink(paste(dir, "/r.sim", sep = ""), append = FALSE)
                cat(paste("r", rnamesX, sep = ""), "\n", sep = "      "); sink() }
   else
     if ("r.sim" %in% FILES) file.remove(paste(dir, "/r.sim", sep = ""))   

   sink(paste(dir, "/otherp.sim", sep = ""), append = FALSE)
   cat("eta", "\n", sep = "      "); sink()

   sink(paste(dir, "/MHinfo.sim", sep = ""), append = FALSE)
   cat("accept.spl.comb", "split", "accept.birth.death", "birth  ", sep = "  ")
   if (nX > 0) cat(paste("beta.block.", 1:nBetaBlocks, sep = ""), "  ", sep = "  ")
   cat("\n", sep = "  ");
   sink()

   if (store$MHb & nrandom > 0){
     sink(paste(dir, "/MHbinfo.sim", sep = ""), append = FALSE)
     cat(paste("b.block.", rep(1:nbBlocks, ncluster), ".", rep(unique.cluster, rep(nbBlocks, ncluster)), "", sep = ""), sep = "  ")
     cat("\n", sep = "  ");
     sink()
   }
   else
     if ("MHbinfo.sim" %in% FILES) file.remove(paste(dir, "/MHbinfo.sim", sep = ""))     

   
   if (store$u){ sink(paste(dir, "/u.sim", sep = ""), append = FALSE)
                 headeru <- paste(" u.", rep(1:prior$kmax, rep(3, prior$kmax)), ".", rep(1:3, prior$kmax), sep = "")
                 headeru[1] <- "    mood"; headeru[2] <- " u.0.0"; headeru[3] <- " u.0.0"
                 cat(headeru, "\n", sep = "     ");
                 sink()
               }
   else
     if ("u.sim" %in% FILES) file.remove(paste(dir, "/u.sim", sep = ""))

   if (store$regresres){sink(paste(dir, "/regresres.sim", sep = ""), append = FALSE)
                        cat(paste("res", rnamesX, sep = ""), "\n", sep = "      "); sink() }
   else
     if ("regresres.sim" %in% FILES) file.remove(paste(dir, "/regresres.sim", sep = ""))      
   
}
