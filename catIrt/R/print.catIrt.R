print.catIrt <-
function(x, digits = max(3, getOption("digits") - 3), ... )
{
   
###############
# Model Stuff #
###############

# If the model is binary, we don't need to let them know about the categories.
  if( x$mod$mod == "brm" )
    cat( "\nModel:\n  ",
         "Binary response model (brm) with ", length(x$cat_theta), " simulees and ", dim(x$full_params)[1],
         " items in the bank.\n\n", sep = "" )

# If the model is polytomous, we should let them know about the categories:
  else
    cat( "\nModel:\n  ",
         "Graded response model (grm) with ", length(x$cat_theta), " simulees, ", dim(x$full_params)[1],
         " items in the bank, and ", dim(x$full_params)[2] - 2, " levels per item.\n\n", sep = "" )
    

#####################
# Termination Stuff #
#####################

  cat("Test Termination:\n")

## FIXED ##
  if( any(x$mod$catTerm$term == "fixed") )
    cat( "  Fixed criterion of ", x$mod$catTerm$n.max, " items.\n", sep = "" )
    
## PRECISION ##
  if( any(x$mod$catTerm$term == "precision") ){
    y <- switch( x$mod$catTerm$p.term$method,
                 threshold = paste("an SEM less than ", x$mod$catTerm$p.term$crit, sep = ""),
                 change    = paste("a change in SEM less than ", x$mod$catTerm$p.term$crit, sep = "") ) 
    cat( "  Variable criterion of ", y, ".\n", sep = "" )
  }
  
## INFO ##
  if( any(x$mod$catTerm$term == "info") ){
    y <- switch( x$mod$catTerm$p.term$method,
                 threshold = paste("Fisher information greater than ", x$mod$catTerm$i.term$crit, sep = ""),
                 change    = paste("a change in Fisher information less than ", x$mod$catTerm$i.term$crit, sep = "") ) 
    cat( "  Variable criterion of ", y, ".\n", sep = "" )
  }
    
## CLASSIFICATION ##
  if( any(x$mod$catTerm$term == "class") ){
    y <- switch( x$mod$catTerm$c.term$method,
                 SPRT = "a Sequential Probability Ratio Test (SPRT)",
                 GLR  = "a Generalized Likelihood Ratio (GLR)",
                 CI   = "a Confidence Interval (CI)" )
    cat( "  Classification using ", y, " criterion.\n", sep = "" )
  } # END if STATEMENT


###################
# Theta Estimates #
###################

# Theta is the most important thing - everything else is in the summary:
  if( length(x$cat_theta) == 1 ){
  	print( summary(x, group = FALSE, ids = 1), digits = digits )
  } else if( length(x$cat_theta) <= 200 ){
    cat( "\nFinal CAT Theta Estimates:\n" )
    print(x$cat_theta, digits = digits)
  } else{
    cat( "\nUse \"catIrt\" element $cat_theta to obtain your Final CAT Theta Estimates.\n")
  } # END ifelse STATEMENT
  
} # END print.catIrt

