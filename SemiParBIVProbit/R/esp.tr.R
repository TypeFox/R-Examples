esp.tr <- function(vrb.st, margin){
 
 
mub <- c("BE") 
 

if( !(margin %in% mub) ){
 
   vrb.st <- ifelse( vrb.st > 28,   28, vrb.st )  
   vrb.st <- ifelse( vrb.st < -17, -17, vrb.st ) 
   vrb    <- exp(vrb.st)
    
}


if( margin %in% mub ){
 
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -16, -16, vrb.st ) 
   vrb    <- plogis(vrb.st)
    
}


    
    
 list(vrb = vrb, vrb.st = vrb.st )  
 
}    