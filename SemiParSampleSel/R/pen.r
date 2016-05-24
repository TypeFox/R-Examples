pen <- function(params, qu.mag, sp, VC, eq1 = "yes"){


if( eq1 == "no" ) l.sp1 <- 0 else l.sp1 <- VC$l.sp1  

    ma1 <- matrix(0,VC$gp1,VC$gp1)   
    if(l.sp1 == 0) EQ1P <- adiag(ma1)
    
    if(l.sp1 != 0){
    ind <- 1:l.sp1
    S1 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S1 <- do.call(adiag, lapply(S1, unlist))
    EQ1P <- adiag(ma1, S1)
                   } 
                     
    ma2 <- matrix(0,VC$gp2,VC$gp2) 
    if(VC$l.sp2 == 0) EQ2P <- adiag(ma2)    
    
    if(VC$l.sp2 != 0){
    ind <- (l.sp1 + 1):(l.sp1 + VC$l.sp2)
    S2 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S2 <- do.call(adiag, lapply(S2, unlist))
    EQ2P <- adiag(ma2, S2)
                   }
                            

if(!is.null(VC$gp3)){

  EQ4P <- EQ5P <- EQ6P <- NULL 
 
    ma3 <- matrix(0,VC$gp3,VC$gp3) 
    if(VC$l.sp3 == 0) EQ3P <- adiag(ma3)    
    
    if(VC$l.sp3 != 0){
    ind <- (l.sp1 + VC$l.sp2 + 1):(l.sp1 + VC$l.sp2 + VC$l.sp3)
    S3 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S3 <- do.call(adiag, lapply(S3, unlist))
    EQ3P <- adiag(ma3, S3)
                   }
        
    
    if(!is.null(VC$gp4)){
    
    ma4 <- matrix(0,VC$gp4,VC$gp4) 
    if(VC$l.sp4 == 0) EQ4P <- adiag(ma4)    
    
    if(VC$l.sp4 != 0){
    ind <- (l.sp1 + VC$l.sp2 + VC$l.sp3 + 1):(l.sp1 + VC$l.sp2 + VC$l.sp3 + VC$l.sp4)
    S4 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
    S4 <- do.call(adiag, lapply(S4, unlist))
    EQ4P <- adiag(ma4, S4)
                   }
    
    
                        }
                        
	    if(!is.null(VC$gp5)){
	    
	    ma5 <- matrix(0,VC$gp5,VC$gp5) 
	    if(VC$l.sp5 == 0) EQ5P <- adiag(ma5)    
	    
	    if(VC$l.sp5 != 0){
	    ind <- (l.sp1 + VC$l.sp2 + VC$l.sp3 + VC$l.sp4 + 1):(l.sp1 + VC$l.sp2 + VC$l.sp3 + VC$l.sp4 + VC$l.sp5)
	    S5 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY=FALSE)
	    S5 <- do.call(adiag, lapply(S5, unlist))
	    EQ5P <- adiag(ma5, S5)
	                     }    
                                }
                                
}else{

    if(VC$margins[2] %in% c("N", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2", "G")) 
      {EQ3P <- 0; EQ4P <- 0;    EQ5P <- NULL  }
    if(VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE" ))                   
      {EQ3P <- 0; EQ4P <- NULL; EQ5P <- NULL  }
    if(VC$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG"))              
      {EQ3P <- 0; EQ4P <- 0;    EQ5P <- 0     }


}



    if( eq1 == "no" ) {
    
 
    if(VC$margins[2] %in% c("N", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2", "G")) EQ4P <- EQ5P <- NULL  
    if(VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE"))                                                             EQ3P <- EQ4P <- EQ5P <- NULL  
    if(VC$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG"))                                         EQ5P <- NULL
    
    
    
    S.h <- adiag(EQ2P, EQ3P, EQ4P, EQ5P)
    
    
    
    
    } else S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P)
    
    
    S.h1 <- 0.5*crossprod(params,S.h)%*%params
    S.h2 <- S.h%*%params
    
    list(S.h = S.h, S.h1 = S.h1, S.h2 = S.h2)
   
         }



















