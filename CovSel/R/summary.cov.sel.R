summary.cov.sel<-function(object, ...){

   at <- attributes(object)
   nn <- at$names      
   alla <- object$covar    
   
   if(length(nn) == 7){
   X.T <- object[[1]]
   Q.0 <- object[[2]]
   Q.1 <- object[[3]]
   metod <- object$method
   
   rem0 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Q.0)){
        if( identical(rem0[j], Q.0[i]) ){
          rem0 <- rem0[-j]
        }
    }
   }
 
   rem1 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Q.1)){
        if( identical(rem1[j], Q.1[i]) ){
          rem1 <- rem1[-j]
        }
    }
   }
  
   cat("\n","Original covariate vector:", "\n",alla,"\n","\n", 
   "Minimal subsets of the covariate vector:","\n","Q.0  = ", Q.0,"\n","Q.1  = ", Q.1,
    "\n","\n", "Removed variables:",
    "\n", "Q.0comp  = ", rem0,"\n", "Q.1comp  = ", rem1, "\n","\n",
    "method  = ",metod,"\n","\n")
    
  }else if(length(nn) == 8){
    if(is.null(object$method)){ X.T <- object[[1]]
                                Q.0 <- object[[2]]
                                Q.1 <- object[[3]]
                               
                                
                                rem0 <- alla
                                
                                for(j in 1:length(alla)){
                                  
                                  for(i in 1:length(Q.0)){
                                    if( identical(rem0[j], Q.0[i]) ){
                                      rem0 <- rem0[-j]
                                    }
                                  }
                                }
                                
                                rem1 <- alla
                                
                                for(j in 1:length(alla)){
                                  
                                  for(i in 1:length(Q.1)){
                                    if( identical(rem1[j], Q.1[i]) ){
                                      rem1 <- rem1[-j]
                                    }
                                  }
                                }
                                
                                cat("\n","Original covariate vector:", "\n",alla,"\n","\n",
                                "Minimal subsets of the covariate vector:","\n",
                                "Q.0  = ", Q.0,"\n","Q.1  = ", Q.1,
                                    "\n","\n", "Removed variables:",
                                    "\n", "Q.0comp  = ", rem0,"\n", "Q.1comp  = ", rem1, "\n","\n",
                                    "regtype = ",object$regtype,"\n","\n")}else{
    
   X.0 <- object[[1]]
   X.1 <- object[[2]]
   Z.0 <- object[[3]]
   Z.1 <- object[[4]]
   metod <- object$method
   
   rem0 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Z.0)){
        if( identical(rem0[j], Z.0[i]) ){
          rem0 <- rem0[-j]
        }
    }
   }
 
   rem1 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Z.1)){
        if( identical(rem1[j], Z.1[i]) ){
          rem1 <- rem1[-j]
        }
    }
   }
 
   cat("\n","Original covariate vector:", "\n",alla,"\n","\n",
   "Minimal subsets of the covariate vector:","\n",
   "Z.0  = ", Z.0,
   "\n","Z.1  = ", Z.1, "\n","\n", "Removed variables:",
    "\n", "Z.0comp  = ", rem0, "\n","Z.1comp  = ", rem1, "\n","\n",
    "method  = ",metod,"\n","\n")} 
}else if(length(nn) == 9){
   
                                      
                                      X.0 <- object[[1]]
                                      X.1 <- object[[2]]
                                      Z.0 <- object[[3]]
                                      Z.1 <- object[[4]]
                                     
                                      
                                      rem0 <- alla
                                      
                                      for(j in 1:length(alla)){
                                        
                                        for(i in 1:length(Z.0)){
                                          if( identical(rem0[j], Z.0[i]) ){
                                            rem0 <- rem0[-j]
                                          }
                                        }
                                      }
                                      
                                      rem1 <- alla
                                      
                                      for(j in 1:length(alla)){
                                        
                                        for(i in 1:length(Z.1)){
                                          if( identical(rem1[j], Z.1[i]) ){
                                            rem1 <- rem1[-j]
                                          }
                                        }
                                      }
                                      
                                      cat("\n","Original covariate vector:", "\n",alla,"\n","\n",
                                      "Minimal subsets of the covariate vector:","\n",
                                      "Z.0  = ", Z.0,
                                          "\n","Z.1  = ", Z.1, "\n","\n", "Removed variables:",
                                          "\n", "Z.0comp  = ", rem0, "\n","Z.1comp  = ", rem1, "\n","\n",
                                          "regtype  = ",object$regtype,"\n","\n")
  

}else if(length(nn) == 13){
   X.T <- object[[1]]
   Q.0 <- object[[2]]
   Q.1 <- object[[3]]
   X.0 <- object[[4]]
   X.1 <- object[[5]]
   Z.0 <- object[[6]]
   Z.1 <- object[[7]]
   metod <- object$method
   
   rem0 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Q.0)){
        if( identical(rem0[j], Q.0[i]) ){
          rem0 <- rem0[-j]
        }
    }
   }
 
   rem1 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Q.1)){
        if( identical(rem1[j], Q.1[i]) ){
          rem1 <- rem1[-j]
        }
    }
   }
   
   rem2 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Z.0)){
        if( identical(rem2[j], Z.0[i]) ){
          rem2 <- rem2[-j]
        }
    }
   }
 
   rem3 <- alla
   
   for(j in 1:length(alla)){
 
    for(i in 1:length(Z.1)){
        if( identical(rem3[j], Z.1[i]) ){
          rem3 <- rem3[-j]
        }
    }
   }
   
   cat("\n","Original covariate vector:", "\n",alla,"\n","\n",
   "Minimal subsets of the covariate vector:","\n",
   "Q.0  = ", Q.0,"\n",
   "Q.1  = ", Q.1,"\n",
   "Z.0  = ", Z.0,
   "\n","Z.1  = ", Z.1, "\n","\n",
   "Removed variables:","\n", "Q.0comp  = ", rem0,"\n", "Q.1comp  = ", rem1,
   "\n", "Z.0comp  = ", rem2, "\n","Z.1comp  = ", rem3, "\n","\n",
   "method  = ",metod,"\n","\n")
   
  }else if(length(nn) == 14){
    X.T <- object[[1]]
    Q.0 <- object[[2]]
    Q.1 <- object[[3]]
    X.0 <- object[[4]]
    X.1 <- object[[5]]
    Z.0 <- object[[6]]
    Z.1 <- object[[7]]
    regtype <- object$regtype
    
    rem0 <- alla
    
    for(j in 1:length(alla)){
      
      for(i in 1:length(Q.0)){
        if( identical(rem0[j], Q.0[i]) ){
          rem0 <- rem0[-j]
        }
      }
    }
    
    rem1 <- alla
    
    for(j in 1:length(alla)){
      
      for(i in 1:length(Q.1)){
        if( identical(rem1[j], Q.1[i]) ){
          rem1 <- rem1[-j]
        }
      }
    }
    
    rem2 <- alla
    
    for(j in 1:length(alla)){
      
      for(i in 1:length(Z.0)){
        if( identical(rem2[j], Z.0[i]) ){
          rem2 <- rem2[-j]
        }
      }
    }
    
    rem3 <- alla
    
    for(j in 1:length(alla)){
      
      for(i in 1:length(Z.1)){
        if( identical(rem3[j], Z.1[i]) ){
          rem3 <- rem3[-j]
        }
      }
    }
    
    cat("\n","Original covariate vector:", "\n",alla,"\n","\n",
    "Minimal subsets of the covariate vector:","\n",
    "Q.0  = ", Q.0,"\n","Q.1  = ", Q.1,"\n",
        "Z.0  = ", Z.0,"\n","Z.1  = ", Z.1, "\n","\n",
        "Removed variables:","\n", "Q.0comp  = ", rem0,"\n", "Q.1comp  = ", rem1,
        "\n", "Z.0comp  = ", rem2, "\n","Z.1comp  = ", rem3, "\n","\n",
        "regtype  = ",regtype,"\n","\n")
  }
}