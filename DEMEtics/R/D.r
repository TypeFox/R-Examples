D<-function(Hs,Ht,values){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           Hs,Ht,values <- D.calc();
#------------------------------------------------------------------------------------------------------------------------------  

          # A function that calculates the D value.
          # The arguments are the Hs-value, the Ht-value and the vector 'values'
          # (a vector that contains the number of individuals that have been 
          # sampled from each population - the sample sizes).
          
          # The Hs- and Ht-values can be calculated by the functions 'Hs()' and 'Ht()'.
           
          # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
          # Molecular Ecology 17,4015-4026.


          n<-length(values)

          ((Ht-Hs)/
                          (1-Hs))*
                                      (n/(n-1))
                                      
                    # Out of Ht, Hs and n, the D value is calculated.                                      
                                                  
}
