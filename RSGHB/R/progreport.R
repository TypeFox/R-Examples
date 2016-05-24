progreport <- function(r, p, a, b, d, f, env)
{

     if (env$gNIV > 0)
     {
          paramRMS <- sqrt(mean(apply(trans(b,env), 1, function(x) x^2)))   
          avgVariance <- mean(apply(trans(b, env), 1, var))
     }
     else
     {
          paramRMS <- NA
          avgVariance <- NA
     }
     
     # Some printing progress to the screen
     if (env$verbose) {
          cat(rep("\n", 128))
          cat("-----------------------------------------------------------\n")
          cat("Iteration: ", r, "\n",sep="\t")
          cat("-----------------------------------------------------------\n")
          # Model Statistics
          mstats <- data.frame(` ` = c("RHO (Fixed):",
                                       "Acceptance Rate (Fixed):",
                                       "RHO (Normal):",
                                       "Acceptance Rate (Normal):",
                                       "Parameter RMS:",
                                       "Avg. Variance:",
                                       "Log-Likelihood:",
                                       "RLH:"),
                               ` ` = as.character(NA), check.names = FALSE, stringsAsFactors = FALSE)
          
          if(env$gFIV > 0) {  
               mstats[1, 2] <- signif(env$rhoF, env$gSIGDIG)
               mstats[2, 2] <- signif(env$acceptanceRateFPerc, env$gSIGDIG)
          }
          
          if(env$gNIV > 0) {
               mstats[3, 2] <- signif(env$rho, env$gSIGDIG)
               mstats[4, 2] <- signif(env$acceptanceRatePerc, env$gSIGDIG)
               mstats[5, 2] <- signif(paramRMS, env$gSIGDIG)
               mstats[6, 2] <- signif(avgVariance, env$gSIGDIG)
          }

          mstats[7, 2] <- signif(sum(log(p)), env$gSIGDIG)
          mstats[8, 2] <- signif(mean(p^(1/env$TIMES)), env$gSIGDIG)
          print(mstats[complete.cases(mstats), , drop = FALSE], right = TRUE, row.names = FALSE)
          
          cat("\n-----------------------------------------------------------\n\n")
          
          # fixed coefficients
          if(env$gFIV > 0)
          {
               print(data.frame(`Fixed Parameters` = paste0(env$gVarNamesFixed, ":"), Estimate = signif(f, env$gSIGDIG), check.names = FALSE), row.names = FALSE)
               cat("\n-----------------------------------------------------------\n\n")
          }
          
          # Normal Coefficients
          if(env$gNIV > 0)
          {
               print(data.frame(`Random Parameters` = paste0(env$gVarNamesNormal, ":"), Estimate = signif(a, env$gSIGDIG), check.names = FALSE), row.names = FALSE)
               cat("\n-----------------------------------------------------------\n")
          }
          
          # Estimated time to completion
          if(r > 1)
          {               
               tpi <- (Sys.time() - env$starttime)/env$gINFOSKIP
               env$starttime <- Sys.time() # makes the forecast based on the most recent iterations
               tleft <- (env$gNCREP + env$gNEREP * env$gNSKIP - r)*tpi
               units(tleft) <- "mins"
               
               cat("Time per iteration:", format(tpi, digits = 3))
               cat("\n")
               cat("Time to completion:", format(tleft, digits = 3))
               cat("\n")
          } else {
               cat("Time per iteration: Calculating...")
               cat("\n")
               cat("Time to completion: Calculating...")
               cat("\n")
          }
          cat("-----------------------------------------------------------\n")
          
          # plotting of the normals  to the screen
          if (env$gNIV > 0 & env$gFIV > 0) {
               
               alphas <- c(r, a, f)
               
          } else if (env$gNIV > 0) {
               
               alphas <- c(r, a)
               
          } else if (env$gFIV > 0) {
               
               alphas <- c(r, f)
               
          }
          
          cr <- rainbow(length(alphas) - 1)
          
          # initial plot
          if (r == 1) {
               xmax <- (env$gNCREP + env$gNEREP * env$gNSKIP) * 1.05
               plot(x = 0, y = 0, main = "Markov Chains", xlim = c(0, xmax), ylim = c(-5, 5), pch = 20, xlab = "Iterations", ylab = "Utility", axes = FALSE, col = "white", cex = 0.5)
               segments(env$gNCREP, -100, env$gNCREP, 100, col = "red", lty = 2, lwd = 2)
               segments(0, 0, env$gNCREP + env$gNEREP * env$gNSKIP, 0, col = "gray", lty = 1, lwd = 1)
               axis(1, at = seq(from = 0, to = env$gNCREP + env$gNEREP * env$gNSKIP, by = floor((env$gNCREP + env$gNEREP * env$gNSKIP) / 10)))
               axis(2, at = -100:100)
          }
          
     
          # Plotting Alphas
          for (i in 2:length(alphas)) points(x = alphas[1], y = alphas[i], pch = 20, col = cr[i - 1], cex = 0.5)
          
          Sys.sleep(0) # forces replotting
     }
     
     # Sink iteration details to the log file.
     detail <- c(r,
                 signif(sum(log(p)), env$gSIGDIG),
                 signif(mean(p^(1/env$TIMES)), env$gSIGDIG),
                 signif(paramRMS, env$gSIGDIG),
                 signif(avgVariance, env$gSIGDIG),
                 if (env$gFIV > 0) {signif(env$acceptanceRateFPerc, env$gSIGDIG)} else {NA},
                 if (env$gNIV > 0) {signif(env$acceptanceRatePerc, env$gSIGDIG)} else {NA})
     
     env$results[["iter.detail"]] <- rbind(env$results[["iter.detail"]], detail)
     
}
