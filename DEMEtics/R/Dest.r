Dest<-function(Hs,Ht,values){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           Hs,Ht,values <- Dest.calc();
#------------------------------------------------------------------------------------------------------------------------------  

          # A function that calculates the Dest value.
          # The arguments are the Hs-value, the Ht-value and the vector 'values'
          # (a vector that contains the number of individuals that have been 
          # sampled from each population - the sample sizes).
          
          # The Hs- and Ht-values can be calculated by the functions 'Hs()' and 'Ht()'.
           
          # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
          # Molecular Ecology 17,4015-4026.


          n<-length(values)

                    # n is the number of populations that have been sampled.

                                      
                    # The harmonic mean is calculated from the sample sizes...
                                                           
          N<-harmonic(values,n)
         
                    # ... and ascribed to the variable N.


                                                       
          Hs.est<-Hs.est(Hs,N)
          
                    # The Hs.est values are calculated and ascribed to the variable 'Hs.est'.
                    # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
                    # Molecular Ecology 17,4015-4026.  

          Ht.est<-Ht.est(Ht,Hs.est,N,n)
          
                    # The Ht.est values are calculated and ascribed to the variable 'Ht.est'.
                    # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
                    # Molecular Ecology 17,4015-4026. 

          ((Ht.est-Hs.est)/
                          (1-Hs.est))*
                                      (n/(n-1))
                                      
                    # Out of Ht.est, Hs.est and n, the Dest value is calculated.                                      
                                                  
}
