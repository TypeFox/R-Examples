trans <- function(b,env)
{
     C <- b
     
     gDIST    <- env$gDIST
     gMAXCOEF <- env$gMAXCOEF
     gMINCOEF <- env$gMINCOEF
     
     for(k in 1:dim(b)[2])
     {
          
          if( gDIST[k] == 2 )	     # positive log-normal
          {
               C[,k] <- exp(b[,k])
          }
          if( gDIST[k] == 3 )	     # negative log-normal
          {
               C[,k] <- -1*exp(b[,k])
          }    
          if( gDIST[k] == 4 )      # censored normal where negative numbers are massed @ 0
          {
               C[,k] <- b[,k]*(b[,k] >= 0)
          }
          if( gDIST[k] == 5 )       # censored normal where positive numbers are massed @ 0
          {
               C[,k] <- b[,k]*(b[,k] <= 0)
          }           
          if( gDIST[k] == 6 )        # Johnson SB
          {	
               C[,k] <- exp(b[,k]) / (1 + exp(b[,k]))
               C[,k] <- (C[,k] * ( gMAXCOEF[k] - gMINCOEF[k] ) ) + gMINCOEF[k]    	
          }             
     }
     
     return(C)
}
