teta.tr <- function(VC, teta.st){
 
 epsilon <- 0.0000001 
 cjg <- c("C0","C180","C90","C270","J0","J180","J90","J270","G0","G180","G90","G270")
 
 if( VC$BivD %in% c("N","AMH","FGM") ) teta.st <- ifelse( abs(teta.st) > 8.75, sign(teta.st)*8.75, teta.st )  
      
 if( VC$BivD %in% cjg ) {
      teta.st <- ifelse( teta.st > 28, 28, teta.st )     # 709, maximum allowed
      teta.st <- ifelse( teta.st < -17, -17, teta.st )   # -20
                         }
      
     
 if(VC$BivD %in% c("N","AMH","FGM") )           teta <- tanh(teta.st)                        
 if(VC$BivD == "F")                             teta <- ifelse( abs(teta.st) < epsilon, epsilon, teta.st )
 if(VC$BivD %in% c("C0", "C180") )              teta <-    exp(teta.st)                   
 if(VC$BivD %in% c("C90","C270") )              teta <- -( exp(teta.st) )                     
 if(VC$BivD %in% c("J0", "J180","G0", "G180") ) teta <-    exp(teta.st) + 1     
 if(VC$BivD %in% c("J90","J270","G90","G270") ) teta <- -( exp(teta.st) + 1 ) 
    
 list(teta = teta, teta.st = teta.st)   
 
   
}    