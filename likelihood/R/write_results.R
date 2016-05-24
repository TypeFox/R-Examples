##################################################
# write_results
# Writes annealing results to a file in a nicely formatted way.
# Parameters:
# results - object returned from anneal
# filename - filename to write to
# print_whole_hist - whether to print the whole likelihood history or an 
# abbreviated version
#
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
#################################################
write_results<-function(results, filename, data=TRUE, print_whole_hist=FALSE) {
  if (!is.character(filename)) stop("filename is not a character string.\n")
   
  out<-file(filename, open="w")

  tryCatch(
  {

  ## Count the number of parameters
  numpars <- 0
  for (i in 1:length(results$par_lo))
    numpars <- numpars + length(results$par_lo[[i]])

  ## Header
  if (results$note=="") {
    cat("Note: [none]\n", sep="", file=out)
  } else {
    cat("Note: ", results$note, "\n", sep="", file=out)
  }
  cat("Total # of Iterations:",results$likeli_hist$iter[[length(results$likeli_hist$iter)]],"\n", sep="\t", file=out)
  cat("Sample size:",length(results$source_data[[1]]),"\n", sep="\t", file=out)
  cat("Maximum likelihood:\t# Parameters:\tAIC corr:\tAIC:\tSlope:\tR2:\n",sep="", file=out)
  cat(results$max_likeli, numpars, results$aic_corr,
      results$aic, results$slope, results$R2, "\n", sep="\t", file=out)
  cat("Model:\t", results$model, "\n", sep="", file=out)
  cat("PDF:\t", results$pdf, "\n", sep="", file=out)
  cat("Support interval range:\t", results$support_interval_range, "\n", sep="", file=out)

  ## Parameter list
  cat("\nVarying parameters:\n", sep="", file = out)
  cat("Label\tMLE\tLower S.I.\tUpper S.I.\tStd Err\tLower Bound\tUpper Bound\tFinal step size\n", sep = "", file=out)
  for (i in 1:length(results$best_pars)) {
    if (length(results$best_pars[[i]]) == 1) {
      cat(names(results$best_pars)[[i]], results$best_pars[[i]],
           results$lower_limits[[i]], results$upper_limits[[i]],
           results$std_errs[[i]],
           results$par_lo[[i]], results$par_hi[[i]], results$par_step[[i]], "\n",
           sep = "\t", file = out)
    }
    else {
      for (j in 1:length(results$best_pars[[i]])) {
        if (length(names(results$best_pars[[i]])) == 0)
          cat(paste(names(results$best_pars)[[i]],j,sep=""),
             results$best_pars[[i]][[j]],
             results$lower_limits[[i]][[j]], results$upper_limits[[i]][[j]],
             results$std_errs[[i]][[j]],
             results$par_lo[[i]][[j]],
             results$par_hi[[i]][[j]],
             results$par_step[[i]][[j]], "\n",
             sep = "\t", file = out)
        else
          cat(names(results$best_pars[[i]])[j],
             results$best_pars[[i]][[j]],
             results$lower_limits[[i]][[j]], results$upper_limits[[i]][[j]],
             results$std_errs[[i]][[j]],
             results$par_lo[[i]][[j]],
             results$par_hi[[i]][[j]],
             results$par_step[[i]][[j]], "\n",
             sep = "\t", file = out)
      }
    }
  }

  cat("\nNon-varying parameters:\n", sep="", file = out)
  cat("Label\tValue\n", sep = "", file=out)
  if (length(var) == 0) cat("[None]", sep = "", file=out)
  else {
    for (i in 1:length(results$var)) {
      if (mode(results$var[[i]]) == "function")
      cat(names(results$var)[[i]], "[function]", "\n", sep = "\t", file = out)
      else if (mode(results$var[[i]]) == "list") {
        for (j in 1:length(results$var[[i]])) {
          cat(paste(names(results$var)[[i]],"$",names(results$var[[i]])[[j]],sep=""),
          results$var[[i]][[j]], "\n",
          sep = "\t", file = out)
        }
      }
      else if (length(results$var[[i]]) == 1)
      cat(names(results$var)[[i]], results$var[[i]], "\n", sep = "\t", file = out)
      else {
        for (j in 1:length(results$var[[i]])) {
          cat(paste(names(results$var)[[i]],j,sep=""),
          results$var[[i]][[j]], "\n",
          sep = "\t", file = out)
        }
      }
    }
  }
  cat("\n", sep="", file = out)
  
  ## Write the initial annealing regime
  cat ("\nAnnealing regime:\n", sep="", file = out)
  cat("Initial temp\tRt\tns\tnt\n", sep="", file = out)
  cat(results$initial_temp, results$temp_red, results$ns, results$nt, sep = "\t", file = out)
  cat("\n", sep="", file = out)
  
  ## Write the variance/covariance matrix, if it was
  ## calculated (anneal ran with hessian = TRUE)
  cat("\nVariance - Covariance Matrix:\n", sep="", file = out)
  if(is.numeric(results$var_covar_mat)) {
    flat_par_names <- NULL

    # Write the parameter names as columns across the top,
    # so watch out for ragged arrays
    for (i in 1:length(results$best_pars)) {
      if (length(results$best_pars[[i]]) == 1) {
        flat_par_names = c(flat_par_names, names(results$best_pars)[i])
      }
      else {
        for (j in 1:length(results$best_pars[[i]]))
          flat_par_names = c(flat_par_names, paste(names(results$best_pars)[i],j,sep=""))
      }
    }
    cat("\t", sep="", file = out)
    cat(flat_par_names, sep="\t", file=out)
    cat("\n", sep="", file = out)
    for (i in 1:length(flat_par_names)) {
      cat(flat_par_names[i], results$var_covar_mat[i,], "\n", sep="\t", file=out)
    }
#    write(results$var_covar_mat, file = out, ncolumns = count, sep = "\t");
  }
  else {
    # Variance / covariance matrix not calculated, so 
    # allow the not calculated message from anneal to be 
    # printed
    cat(results$var_covar_mat, file = out, sep = "")
    cat("\n", sep="", file = out)
  }


  ## Source data - I'm not using write.table because it gives a
  ## warning message about appending a table to a file, and that's
  ## just stupid and annoying.
  cat("\nSource data\n", sep="", file = out)
  if (data) {
    for (i in 1:length(results$source_data))
      cat(names(results$source_data)[[i]], "\t", sep="", file = out)
    cat("\n", sep="", file=out)
    for (i in 1:length(results$source_data[[1]])) {
      for (j in 1:length(results$source_data)) {
        cat(results$source_data[[j]][[i]], "\t", sep="", file = out)
      }
      cat("\n", sep="", file=out)
    }
  } else {
    cat("[Ommitted from results file]\n", sep="", file = out)
  }

  ## Likelihood history
  cat("\nLikelihood History\n", sep = "", file = out)
  if (print_whole_hist) {
    cat("Iteration","Temperature","Likelihood",
         colnames(results$likeli_hist)[4:length(colnames(results$likeli_hist))],
         "\n", sep = "\t", file = out)
    for (i in 1:length(results$likeli_hist$temp)) {
      cat(results$likeli_hist$iter[i],
        results$likeli_hist$temp[i],
        results$likeli_hist$likeli[i],sep = "\t", file = out)
      for (j in 4:dim(results$likeli_hist)[2])
        cat("\t",results$likeli_hist[[j]][i],sep = "", file = out)
      cat("\n", sep="", file=out)
    }
  } else {   
    cat("Iteration\tTemperature\tLikelihood\n", sep = "", file = out)
    for (i in 1:length(results$likeli_hist$temp)) {
      # Cut down on the number of results stored - print every thousandth
      # instead of every hundredth, along with all changes in likelihood
      if (!(results$likeli_hist$iter[i] %% 100 == 0 &&
          results$likeli_hist$iter[i] %% 1000 != 0))
      cat(results$likeli_hist$iter[i],
          results$likeli_hist$temp[i],
          results$likeli_hist$likeli[i],  "\n", sep = "\t", file = out)
    }
  }

  }, finally=close(out))
}