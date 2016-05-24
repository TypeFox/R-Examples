Hardy.Weinberg <- function(tab2,l){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#                 tab2 <- Bootstrapping.Chao(), Bootstrapping(), Bootstrapping.Gst();

# Output:
#                 HWE -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------  

          # This function calculates if all populations are in HWE for the
          # actual locus l.

tab2.pop <- split(tab2,tab2$population)

          # The data are splitted so that those belonging to different populations
          # are separated.
          
number.pops <- length(tab2.pop)

          # The number of populations for which data for the actual locus were
          # obtained.
          

Hardy1 <- lapply(tab2.pop,Hardy.calc)
Hardy <- do.call(c,Hardy1)

          # This vector is filled with the p.values for the several
          # populations that give the probabilities that the populations are
          # in HWE for the actual locus.


                          
# The p.values (probabilities if HWE=TRUE) are combined. These are
# obtained from a 10.000 fold repeated Monte Carlo simulation.

HWE <- ifelse(all(Hardy>0.05),TRUE,FALSE)        

          # If all populations are in HWE for this locus, HWE is set as TRUE,
          # otherwise it is set as FALSE. 
          
assign("HWE",HWE,pos = DEMEtics.env)          

          # This end result is assigned to the workspace.
          
}
