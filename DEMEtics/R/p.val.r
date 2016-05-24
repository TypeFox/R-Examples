p.val <- function(empirical.value, bootstrapped.values){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#          empirical.value,bootstrapped.values <- all.pops.Dest.Chao(), all.pops.Dest(), all.pops.Gst(), pair.pops.Gst(),
#                                                 pair.pops.Dest.Chao(), pair.pops.Dest();

# Output:
#          p.value -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------  

          # Function that enables to assign a p-value to an empirical value
          # when a bootstrap procedure has been carried out.
          # This procedure is limited to one.sided tests, when the alternative
          # hypothesis is 'larger than'.

          bt <- length(bootstrapped.values)

          # bt is the number of times the values were bootstrapped.

    p.value <- (1+sum(bootstrapped.values >= empirical.value))/(bt+1)
    
          # The p.values are calculated accoring to 
          # Manly BFJ. (1997). Randomization, bootstrap and Monte Carlo methods 
          # in biology. (Chapman & Hall, London [u.a.]), p. 62.
          
          # bt is the number of repetitions in the bootstrapping process.
          
    assign("p.value",p.value,pos = DEMEtics.env)
    
          # The p.value is assigned to the workspace so that it can be obtained
          # for further calculations.
          
}
    

    
