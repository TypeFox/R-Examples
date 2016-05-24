print.summary.catIrt <-
function(x, digits = max(3, getOption("digits") - 3), ... ){
  
# Note: Prints "Group Level Summary" and "Individual Level Summary"
#       only if both are select.  Otherwise, it prints neither.

#~~~~~~~~~~~~~~~~~~~~~#
# Group Level Summary #
#~~~~~~~~~~~~~~~~~~~~~#
  y <- x$cat_run

  if( x$group == TRUE ){
	
    if( !( is.null(x$ids) | any(x$ids == "none" ) ) ){
      cat("\n   *************************\n")
      cat(  "   ** Group Level Summary **\n")
      cat(  "   *************************\n")
    } # END if STATEMENT
  
###############
# Model Stuff #
###############

# If the model is binary, we don't need to let them know about the categories:
    if( y$mod$mod == "brm" )
      cat( "\nModel:\n  ",
           "Binary response model (brm) with ", x$theta$d_stat["n"], " simulees and ", x$items$bank["n.items"],
           " items in the bank.\n\n\n", sep = "" )

# If the model is polytomous, we should let them know about the categories:
    else
      cat( "\nModel:\n  ",
           "Graded response model (grm) with ", x$theta$d_stat["n"], " simulees, ", x$items$bank["n.items"],
           " items in the bank, and ",x$items$bank["n.params"], " levels per item.\n\n\n", sep = "" )

###############
# Theta Stuff #
###############
    cat("## Summary of Theta: ##\n\n")
      
    if( any(rownames(x$theta$sum) == "true") ){
      cat( "True Theta:\n", sep = "" )
      cat( "    Mean = ", m <- round( x$theta$sum["true", "mean"], digits ), "   ", sep = "")
      cat(     "SD  = ",       round( x$theta$sum["true", "sd"],   digits ), "\n", sep = "")
      cat( "    Min  = ", n <- round( x$theta$sum["true", "min"],  digits ),
           rep(" ", nchar(m) - nchar(n) + 3), sep = "")
      cat(     "Max = ",       round( x$theta$sum["true", "max"],  digits ), "\n",  sep = "")
    } # END if STATEMENT    
  
    cat( "Full Bank Theta (using ", y$mod$catTerm$score, " Estimation):\n",  sep = "")
    cat( "    Mean = ", m <- round( x$theta$sum["tot", "mean"], digits ), "   ",    sep = "")
    cat(     "SD  = ",       round( x$theta$sum["tot", "sd"],   digits ), "\n",    sep = "")
    cat( "    Min  = ", n <- round( x$theta$sum["tot", "min"],  digits ),
         rep(" ", nchar(m) - nchar(n) + 3),    sep = "")
    cat(     "Max = ",       round( x$theta$sum["tot", "max"],  digits ), "\n",     sep = "")

    cat( "CAT Theta (using ", y$mod$catMiddle$score, " Estimation):\n",    sep = "")
    cat( "    Mean = ", m <- round( x$thet$sum["cat", "mean"], digits ), "   ", sep = "")
    cat(     "SD  = ",       round( x$thet$sum["cat", "sd"],   digits ), "\n", sep = "")
    cat( "    Min  = ", n <- round( x$thet$sum["cat", "min"],  digits ),
         rep(" ", nchar(m) - nchar(n) + 3), sep = "")
    cat(     "Max = ",       round( x$thet$sum["cat", "max"],  digits ), "\n",  sep = "")

    if( any(rownames(x$theta$sum) == "true") ){
      cat("CAT Theta vs. True Theta:\n")
      cat( "    Rsq = ",                round(x$theta$e_stat["rsq"],  digits), "\n",  sep = "")
      cat( "    Avg Sign Err = ",  m <- round(x$theta$e_stat["ase"],  digits), "   ",  sep = "")
      cat(     "MSE  = ",               round(x$theta$e_stat["mse"],  digits), "\n",  sep = "")
      cat( "    Avg Abs Err  = ",  n <- round(x$theta$e_stat["aae"],  digits),
           rep(" ", nchar(m) - nchar(n) + 3),  sep = "")
      cat(     "RMSE = ",               round(x$theta$e_stat["rmse"], digits), "\n", sep = "")
    } # END if STATEMENT    
      
    cat("CAT Theta vs. Full Bank Theta:\n")
    cat( "    Rsq = ",                round(x$theta$d_stat["rsq"],  digits), "\n",  sep = "")
    cat( "    Avg Sign Diff = ", m <- round(x$theta$d_stat["asd"],  digits), "   ",  sep = "")
    cat(     "MSD  = ",               round(x$theta$d_stat["msd"],  digits), "\n",  sep = "")
    cat( "    Avg Abs Diff  = ", n <- round(x$theta$d_stat["aad"],  digits),
         rep(" ", nchar(m) - nchar(n) + 3),  sep = "")
    cat(     "RMSD = ",               round(x$theta$d_stat["rmsd"], digits), "\n\n\n", sep = "")
  

#############
# SEM Stuff #
#############
    cat("## Summary of SEM: ##\n\n")
      
    cat( "SEM of Full Bank Theta:\n" )
    cat( "    Mean = ", m <- round( x$sem$sum["tot", "mean"], digits ), "   ", sep = "" )
    cat(     "SD  = ",       round( x$sem$sum["tot", "sd"],   digits ), "\n", sep = "")
    cat( "    Min  = ", n <- round( x$sem$sum["tot", "min"],  digits ),
         rep(" ", nchar(m) - nchar(n) + 3), sep = "")
    cat(     "Max = ",       round( x$sem$sum["tot", "max"],  digits ), "\n",  sep = "" )
  
    cat( "SEM of CAT Theta:\n", sep = "")
    cat( "    Mean = ", m <- round( x$sem$sum["cat", "mean"], digits ), "   ", sep = "")
    cat(     "SD  = ",       round( x$sem$sum["cat", "sd"],   digits ), "\n", sep = "" )
    cat( "    Min  = ", n <- round( x$sem$sum["cat", "min"],  digits ),
         rep(" ", nchar(m) - nchar(n) + 3), sep = "")
    cat(     "Max = ",       round( x$sem$sum["cat", "max"],  digits ), "\n",  sep = "" )
   
    cat("CAT SEM vs. Full Bank SEM:\n")
    cat( "    Rsq = ",                round( x$sem$d_stat["rsq"], digits ), "\n",   sep = "" )
    cat( "    Avg Sign Diff = ", m <- round( x$sem$d_stat["asd"], digits ), "   ",   sep = "" )
    cat(     "MSD  = ",               round( x$sem$d_stat["msd"], digits ), "\n",   sep = "" )
    cat( "    Avg Abs Diff  = ", n <- round( x$sem$d_stat["aad"], digits ),
         rep(" ", nchar(m) - nchar(n) + 3), sep = "")
    cat(     "RMSD = ",               round( x$sem$d_stat["rmsd"], digits ), "\n\n\n", sep = "" )

  
###############
# Items Stuff #
###############
    cat("Number of Items Administered in the CAT:\n")
    cat( "    Mean = ", round( x$items$sum["mean"], digits ), "\n", sep = "" )
    cat( "    SD   = ", round( x$items$sum["sd"],   digits ), "\n", sep = "" )
    cat( "    Min  = ", round( x$items$sum["min"],  digits ), "\n", sep = "" )
    cat( "    Max  = ", round( x$items$sum["max"],  digits ), "\n", sep = "" )


########################
# Classification Stuff #
########################
    if( any( y$mod$catTerm$term == "class" ) ){
      cat("\nNumber of Simulees in Each Category:")
      print( x$freq$categ )
    }
  

#####################
# Termination Stuff #
#####################
    cat( "\nNumber of Simulees with Each Termination Criterion:" )
    print( x$freq$term )
  

##########################
# Total Bank Information #
##########################

    cat("\nMax 'Expected Bank Information' is ", round( unique(max(x$bank_info$info)), digits ),
        " at a theta of ", round( x$bank_info$eval[ which.max(x$bank_info$info) ], digits ),
        " with an associated SEM of ", round( x$bank_info$sem[ which.max(x$bank_info$info) ], digits ),
        ".\n\n", sep = "" )
        
   } # END if STATEMENT   


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Individual Level Summaries #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# If you want individual summaries - print individual summaries:
  if( !( is.null(x$ids) | any(x$ids == "none" ) ) ){
      
    if( x$group == TRUE ){    
      cat("\n")
      cat("   ********************************\n")
      cat("   ** Individual Level Summaries **\n")
      cat("   ********************************\n")
    } # END if STATEMENT
    
    resp.all <- unlist( lapply(y$cat_indiv, FUN = function(z) z$cat_resp) )
    resp.min <- min(resp.all)
    resp.max <- max(resp.all)
  
# For all people do the following:
# - pull out the id number,
# - give a header of the total bank and CAT bank ability estimates,
# - let people know how the selection criterion, ability estimator, and termination criterion,
# - provide an indication as to the items, responses, and CAT theta updates.
    for( i in seq_along(x$ids) ){
    
################
# Header Stuff #
################

## Summary of Simulee Statement ##
      if( length(x$ids) == 1 ){
      	id.i <- x$ids[i]
        cat("\n## Summary of Simulee: ##\n\n", sep = "")
      } else{
        id.i <- x$ids[i]
        cat("\n## Summary of Simulee ", id.i, ": ##\n\n", sep = "")
      } # END ifelse STATEMENT


## Printing of Theta and Theta Estimates ##
      if( length(y$true_theta) >= 1 & all( !is.na(y$true_theta) ) )
        cat("True Theta: ",      round(y$true_theta[id.i], digits), "\n",  sep = "")
        
      cat(  "Full Bank Theta: ", round(y$tot_theta[id.i], digits), "\n",   sep = "")
      cat(  "Final CAT Theta: ", round(y$cat_theta[id.i], digits), "\n\n", sep = "")
      
##########################
# Generalizing CAT Stuff #
##########################

############
# catStart #
############
      start.i   <- y$mod$catStart
      middle.i  <- y$mod$catMiddle
      resp.i    <- y$cat_indiv[[id.i]]$cat_resp
      
      a <- start.i$n.start
      b <- which( { c(FALSE, diff(resp.i) != 0) |
      	            ( (resp.i != resp.min) & (resp.i != resp.max) ) } )[1]
      	            
# Note: b finds the first instance that:
#       - the response changes OR
#       - the response isn't equal to the min/max possible response.

# a) n.start in catStart
# b) the first change in the response (min -- not min or max -- not max)
# IF: i)    MLE out and score == "MLE" --> then b regardless.
#     ii)   MLE out and score != "MLE" --> then minimum of a or b.
#     iii) !MLE out and score == "MLE" --> then maximum of a or b.
#     iv)  !MLE out and score != "MLE" --> then a regardless.
      
      s.it     <- { if( (start.i$leave.after.MLE)  & (middle.i$score == "MLE") ){
      	              b
      	            } else if( (start.i$leave.after.MLE)  & (middle.i$score != "MLE") ){
      	              min(a, b)
      	            } else if( !(start.i$leave.after.MLE) & (middle.i$score == "MLE") ){
      	              max(a, b)
      	            } else{
      	              a
      	            } # END ifelse STATEMENTS
      	          }

      s.select <- switch( y$mod$catStart$select,
                          "UW-FI"  = "by maximizing Fisher information at",
                          "LW-FI"  = "by maximizing Likelihood-Weighted Fisher information",
                          "PW-FI"  = "by maximizing Posterior-Weighted Fisher information",
                          "FP-KL"  = "by maximizing Fixed Pointwise Kullback-Leibler divergence at",
                          "VP-KL"  = "by maximizing Variable Pointwise Kullback-Leibler divergence at",
                          "FI-KL"  = "by maximizing Fixed Integrated Kullback-Leibler divergence at",
                          "VI-KL"  = "by maximizing Variable Integrated Kullback-Leibler divergence at",
                          "random" = "randomly" )
                       
      if( y$mod$catStart$select != "random" ){
        s.at   <- switch( y$mod$catStart$at,
                          theta  = " theta",
                          bounds = " (a) fixed point(s)" )
      } else{
        s.at   <- NULL
      } # END ifelse STATEMENT
                       
      s.score  <- switch( y$mod$catStart$score,
                          fixed  = paste("a fixed value of ", y$mod$catStart$init.theta, sep = ""),
                          step   = paste("increasing or decreasing theta by ", y$mod$catStart$step.size, sep = ""),
                          random = "randomly",
                          WLE    = "Weighted Likelihood (WL) estimation",
                          BME    = "Bayesian Modal (BM) estimation",
                          EAP    = "Expected A Posteriori (EAP) estimation" )

#############
# catMiddle #
#############
        
      m.select <- switch( y$mod$catMiddle$select,
                          "UW-FI"  = "by maximizing Fisher information at",
                          "LW-FI"  = "by maximizing Likelihood-Weighted Fisher information",
                          "PW-FI"  = "by maximizing Posterior-Weighted Fisher information",
                          "FP-KL"  = "by maximizing Fixed Pointwise Kullback-Leibler divergence at",
                          "VP-KL"  = "by maximizing Variable Pointwise Kullback-Leibler divergence at",
                          "FI-KL"  = "by maximizing Fixed Integrated Kullback-Leibler divergence at",
                          "VI-KL"  = "by maximizing Variable Integrated Kullback-Leibler divergence at",
                          "random" = "randomly" )
                       
      if( y$mod$catMiddle$select != "random" ){
        m.at   <- switch( y$mod$catMiddle$at,
                          theta  = " theta",
                          bounds = " (a) fixed point(s)" )
      } else{
        m.at   <- NULL
      } # END ifelse STATEMENT
                       
      m.score  <- switch( y$mod$catMiddle$score,
                          MLE = "Maximum Likelihood (ML) estimation",
                          WLE = "Weighted Likelihood (WL) estimation",
                          BME = "Bayesian Modal (BM) estimation",
                          EAP = "Expected A Posteriori (EAP) estimation" )

###########
# catTerm #
###########
                       
# For EAP, we terminate slightly differently:
      if( y$mod$catMiddle$score == "EAP" ){
        sem <- "posterior standard deviation"
      } else{
        sem <- "standard error of measurement"
      } # END ifelse STATEMENT
                       
      term    <- switch( y$cat_term[id.i],
                         total       = "all of the items in the bank.",
                         fixed       = paste("a fixed number of ", y$mod$catTerm$n.max, " items.", sep = ""),
                         se_thresh   = paste("the ", sem, " was below ", y$mod$catTerm$p.term$crit, ".", sep = ""),
                         se_change   = paste("the ", sem, " changed less than ", y$mod$catTerm$p.term$crit, " between two successive items.", sep = ""),
                         fi_thresh   = paste("Fisher Information was above ", y$mod$catTerm$i.term$crit, ".", sep = ""),
                         fi_change   = paste("Fisher Information changed less than ", y$mod$catTerm$i.term$crit, " between two successive items.", sep = ""),
                         class       = paste("a classification decision of '", y$cat_categ[id.i], "' was made.", sep = "") )

######################
# Printing CAT Stuff #
######################
      cat("At the beginning of the CAT (first ", s.it, " items):", sep = "")
      cat("\n  Items were selected ", s.select, s.at, ".", sep = "")
      cat("\n  Ability was estimated after each step by ", s.score, ".\n\n", sep = "")
      
      cat("For the rest of the CAT:")
      cat("\n  Items were selected ", m.select, m.at, ".", sep = "")
      cat("\n  Ability was estimated after each step by ", m.score, ".\n\n", sep = "")
      
      cat("The CAT was terminated after ", term, "\n", sep = "")

########################
# Specific Items Stuff #
########################

      indiv.i <- y$cat_indiv[[id.i]]
     

      mat1 <- c(NA, 1:nrow(indiv.i$cat_params))                    # the CAT items
      mat2 <- rbind( c(NA, indiv.i$cat_params[ , 1]),              # the bank items
                     c(NA, indiv.i$cat_resp) )                     # the responses
      mat3 <- rbind( round(indiv.i$cat_theta, 2),                  # the ability estimates
                     as.character(round(indiv.i$cat_sem, digits))  # the SEM
                    )
      
# I made sure that we "as.charactered" it to prevent the item numbers/responses from having lots of 0s
      mat <- rbind(mat1, mat2, mat3)
      
      rownames(mat) <- c("CAT Number", "Bank Number", "Response", "Estimate", "SEM")
      colnames(mat) <- rep("", ncol(mat))
      
      print(format(mat, justify = "right"), quote = FALSE)
      cat("\n\n")
      
     } # END for LOOP
   } # END if STATEMENT
      
} # END print.summary.catIrt FUNCTION