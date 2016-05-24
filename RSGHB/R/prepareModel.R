prepareModel <- function(env)
{
     
     # Global variable setting
     env$gNP        <- length(unique(env$choicedata$ID))       # Number of people
     
     env$Choice     <- env$choicedata$Choice
     
     env$gNOBS      <- dim(env$choicedata)[1]
     
     env$TIMES      <- matrix(0,nrow=env$gNP,ncol=1)     		# Number of observations for each person
     env$TIMES[,1]  <- aggregate(env$choicedata$ID,by=list(env$choicedata$ID),length)[,2]
     
     env$gIDS       <- unlist(as.vector(mapply(rep,1:env$gNP,env$TIMES)))     # index map for individual to observation
     
     env$respIDs    <- unique(env$choicedata$ID)  
     
     # Matrix initialization
     # A must have NIV columns and 1 row
     # B must have NP columns and NIV rows
     # Dmat must have NIV columns and NIV rows and be symmetric
     
     if(length(env$gVarNamesNormal) > 0)
     {
          env$A <- matrix(0,nrow=env$gNIV,ncol=1)
          env$B <- matrix(0,nrow=env$gNP,ncol=env$gNIV)
          env$Dmat <- env$priorVariance * diag(env$gNIV)
 
          env$A[,1] <- env$svN        # initialize to analyst specified starting values
     
          env$B <- 1 + env$B
          env$B <- env$B * matrix(t(env$A),nrow=env$gNP,ncol=env$gNIV,byrow=T)
     }
     
}
