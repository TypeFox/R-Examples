summary.MaxEntICA.BinBin <- function(object, ..., Object){
 
  options(digits = 4)
  
  if (missing(Object)){Object <- object} 
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# Maximum entropy distribution of vector of potential outcomes")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  print(Object$Vector_p)

  cat("\n\n# H_max (entropy of p*) ")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat(Object$H_max, "\n")
  
    
  cat("\n\n# R2_H results")
  cat("\n#~~~~~~~~~~~~~\n\n")  
  cat(Object$R2_H, "\n\n")
  }
