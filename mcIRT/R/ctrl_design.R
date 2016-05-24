ctrl_design <-
function(design,aDD,gr,TYPE=TYPE)  
{

  
  if(TYPE=="NRM")  
  {
    # a dif design for each parameter group is committed
    if(length(design) != 2){stop("The design list must contain 2 matrices! The first for the zeta parameters, the second for the lambda parameters.")}
  } else if(TYPE=="NLM") {
                          # a dif design for each parameter group is committed
                          if(length(design) != 4){stop("The design list must contain 4 matrices!   alpha, beta, zetas, lambdas!")}

                         }
  

    # groups
    ngru <- sapply(design,function(x) nrow(x) != nlevels(gr))
    if(any(ngru)){stop("The number of rows must equal the number of groups!")}
      
    # items
    nite <- sapply(design,function(x) ncol(x) != length(aDD))
    if(any(ngru)){stop("The number of colums must equal the number of items!")}
    
    # numbers inside the matrices
    ctrl4 <- sapply(design,function(toob)
    {
      any(toob > nlevels(gr)) | any(toob < 1)
      
    })
    
    if(any(ctrl4)){stop("The numbers inside the matrix should refer to the groups!")}
    
    
    # numbers inside the matrices2
    
    design[[1]][1,] <- 1 
    design[[2]][1,] <- 1
    
    
    ctrl5 <- sapply(design,function(toob)
    {
      ee1 <- sapply(1:nlevels(gr),function(wron)
      {
        any(toob[wron,] > wron)
        
      })
      ee1
      
    })
    if(any(ctrl5)){stop("Check your matrix! At least one line refers to a group with higher number than the line itself!")}

  

## in this section the design is restructured to fit into the grDM function, and to avoid inconsistencies
## for example, it is not possible, that in the 3rd line of the matrix of item 1 (e.g.) a 2 is listed while in the first and second line a 1 is listed. so it is assumed that one means, that the parameter estimation for the third group should be the same as in the second. while the estimation for the second group is like the one in the first, the third MUST also be the same. so this means that the 2 will be replaced by an 1.

  
design_neu <- lapply(design,function(xy)
                {
                
                  eindesign <- sapply(1:nlevels(gr),function(i)
                            {
                             if(i == 1){
                               xy[i,] 
                               } else {
                                      bz <- xy[i,]
                                      
                                      zeilE<- mapply(function(aa,bb)
                                                {
                                                xy[aa,bb]
                                                },aa=bz,bb=1:length(bz))
                                      
                               
                                      }
                              
                            #
                            })

                t(eindesign)
  
                })

design_neu  
}
