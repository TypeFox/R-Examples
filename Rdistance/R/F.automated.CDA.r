F.automated.CDA <- function (detection.data, transect.data, 
                             w.lo=0, w.hi=max(dist),
                             likelihoods=c("halfnorm", "hazrate", "uniform", "negexp", "Gamma"),
                             series=c("cosine", "hermite", "simple"), expansions=0:3, warn=TRUE,
                             area=1, ci=0.95, R=500, by.id=FALSE, plot.bs=FALSE,                          
                             plot=TRUE, ...){
  
  
  
  # Stop and print error if a 'dist' column isn't provided in detection.data
  if(!("dist" %in% names(detection.data)))
    stop("There is no column named 'dist' in your detection.data.")
  
  # extract distance vector from distdata
  dist <- detection.data$dist
  
  
  # function to save results
  f.save.result <- function(results, dfunc, like, ser, expan, 
                            plot) {
    esw <- ESW(dfunc)
    if (!is.na(esw) & (esw > dfunc$w.hi)) {
      scl.ok <- "Not ok"
      scl.ok.flag <- 0
      aic <- NA
    }
    else {
      scl.ok <- "Ok"
      scl.ok.flag <- 1
      aic = AIC(dfunc)
    }
    conv <- dfunc$convergence
    if (conv != 0) {
      if (conv == -1) {
        conv.str <- "Bad"
      }
      else {
        conv.str <- "No"
      }
      aic <- NA
      scl.ok <- "NA"
      scl.ok.flag <- NA
    }
    else {
      conv.str <- "Yes"
    }
    results <- rbind(results, data.frame(like = like, series = ser, 
                                         expansions = expan, converge = conv, scale = scl.ok.flag, 
                                         aic = aic))
    if (nchar(like) < 8) 
      sep1 <- "\t\t"
    else sep1 <- "\t"
    cat(paste(like, sep1, ser, "\t", expan, "\t", conv.str, 
              "\t\t", scl.ok, "\t", round(aic, 4), sep = ""))
    if (plot) {
      plot(dfunc)
      k <- readline(" Next?[entr=y,n]")
      if (length(k) == 0) 
        k <- "y"
    }
    else {
      cat("\n")
      k <- "y"
    }
    list(results = results, k = k)
  }
  
  
  # Fit detection functions (F.dfunc.estim appears 4 times below)
  
  wwarn <- options()$warn
  options(warn = -1)
  fit.table <- NULL
  cat("Likelihood\tSeries\tExpans\tConverged?\tScale?\tAIC\n")
  for (like in likelihoods) {
    if (like == "Gamma") {
      dfunc <- F.dfunc.estim(dist, likelihood = like, w.lo = w.lo, 
                             w.hi = w.hi, ...)
      ser <- ""
      expan <- 0
      fit.table <- f.save.result(fit.table, dfunc, like, 
                                 ser, expan, plot)
      cont <- fit.table$k
      fit.table <- fit.table$results
    }
    else {
      for (expan in expansions) {
        if (expan == 0) {
          ser <- "cosine"
          dfunc <- F.dfunc.estim(dist, likelihood = like, 
                                 w.lo = w.lo, w.hi = w.hi, expansions = expan, 
                                 series = ser, ...)
          fit.table <- f.save.result(fit.table, dfunc, 
                                     like, ser, expan, plot)
          cont <- fit.table$k
          fit.table <- fit.table$results
        }
        else {
          for (ser in series) {
            dfunc <- F.dfunc.estim(dist, likelihood = like, 
                                   w.lo = w.lo, w.hi = w.hi, expansions = expan, 
                                   series = ser, ...)
            fit.table <- f.save.result(fit.table, dfunc, 
                                       like, ser, expan, plot)
            cont <- fit.table$k
            fit.table <- fit.table$results
            if (cont == "n") 
              break
          }
        }
        if (cont == "n") 
          break
      }
      if (cont == "n") 
        break
    }
  }
  if (sum(fit.table$converge != 0) > 0) 
    cat("Note: Some models did not converge or had parameters at their boundaries.\n")
  fit.table$aic <- ifelse(fit.table$converge == 0, fit.table$aic, 
                          Inf)
  fit.table <- fit.table[order(fit.table$aic), ]

  dfunc <- F.dfunc.estim(dist, likelihood = fit.table$like[1], 
                         w.lo = w.lo, w.hi = w.hi, expansions = fit.table$expansions[1], 
                         series = fit.table$series[1], ...)
  if (plot) {
    plot(dfunc)
    mtext("BEST FITTING FUNCTION", side = 3, cex = 1.5, line = 3)
  }

  
  abund <- F.abund.estim(dfunc, detection.data=detection.data, transect.data=transect.data,
                         area=area, ci=ci, R=R, plot.bs=plot.bs, by.id=by.id)

#   }
  cat("\n\n---------------- Final Automated CDS Abundance Estimate -------------------------------\n")
  print(abund)
  options(warn = wwarn)
  abund
}