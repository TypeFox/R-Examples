test.Null.warning <-
function(x, a1a3.pval){

  if(!x & is.null(a1a3.pval)){
   cat("\n")
   cat("  ##############################################\n")
   cat("  ###                WARNING!                ###\n")
   cat("  ##############################################\n\n")   


   cat("  ### Not enough evidence to reject the      ###\n")
   cat("  ### hypothesis test of:                    ###\n") 
   cat("  ###                                        ###\n")        
   cat("  ### H_0 : No marker-by-treatment           ###\n")
   cat("  ###              interaction               ###\n")
   cat("  ###                                        ###\n") 
   cat("  ### Inference for Theta may be unreliable! ###\n\n")
   cat("  ##############################################\n\n") 

  }else if(!x & !is.null(a1a3.pval)){

   cat("\n")
   cat("  ##############################################\n")
   cat("  ###                WARNING!                ###\n")
   cat("  ##############################################\n\n")   


   cat("  ### Not enough evidence to reject both     ###\n")
   cat("  ### hypothesis tests of:                   ###\n") 
   cat("  ### H_0 :                                  ###\n")        
   cat("  ### 1. No marker-by-treatment interaction  ###\n")
   cat("  ###                                        ###\n")  
   cat("  ### 2. The marker positivity threshold     ###\n") 
   cat("  ###    is outside marker bounds            ###\n")  
   cat("  ###                                        ###\n")  
   cat("  ### Inference for Theta may be unreliable! ###\n\n")
   cat("  ##############################################\n\n") 
  }
}
