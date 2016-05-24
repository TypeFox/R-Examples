eta.tr <- function(vrb.st, margin){

mupos <- c("LN","WEI","iG","GA","DAGUM","SM","FISK")
mub   <- c("BE")

if(margin %in% mupos){
   vrb.st <- ifelse( vrb.st > 28,   28, vrb.st )  
   vrb.st <- ifelse( vrb.st < -17, -17, vrb.st ) 
}

if(margin %in% mub){
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -16, -16, vrb.st ) 
}


vrb.st  
 
}    