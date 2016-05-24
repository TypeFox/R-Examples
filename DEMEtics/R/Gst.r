Gst<-function(Hs,Ht,values){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           Hs,Ht,values <- Gst.calc();
#------------------------------------------------------------------------------------------------------------------------------  

          # A function that calculates the Gst value.
          # The arguments are the Hs-value, the Ht-value and the vector 'values'
          # (a vector that contains the number of individuals that have been 
          # sampled from each population - the sample sizes).
          
          # The Hs- and Ht-values can be calculated by the functions 'Hs()' and 'Ht()'.
           
          # See: Jost L. (2008). Gst and its relatives do not measure differentiation. 
          # Molecular Ecology 17,4015-4026.


          (Ht-Hs)/Ht
          
                    # Out of Ht.est, Hs.est and n, the Gst value is calculated.                                      
                       
          }

