## ------------
## Union Class
## ------------

setClassUnion("covKernel", 
              c("covTensorProduct", "covIso", 
                "covAffineScaling", "covScaling", 
                "covUser")) #, "covAdditive0"))
