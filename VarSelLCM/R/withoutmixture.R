setGeneric ( name= "withoutmixture",  def = function(obj){ standardGeneric("withoutmixture")})
## Pour les variables continues
setMethod( f = "withoutmixture", 
           signature(obj="VSLCMresultsContinuous"), 
           definition = function(obj){
             obj@model@omega <-  as.numeric(obj@model@omega)
             names(obj@model@omega) <- colnames(obj@data@data)
             # On remet les valeurs manquantes
             for (j in 1:obj@data@d){
               if (any(obj@data@notNA[,j]==0))
                 obj@data@data[which(obj@data@notNA[,j]==0),j] <- NA
             }
             obj@param@pi <- 1
             obj@param@mu <- matrix(NA, obj@data@d, 1)
             obj@param@sd <- matrix(NA, obj@data@d, 1)
             rownames(obj@param@mu)  <-   colnames(obj@data@data)
             rownames(obj@param@sd)  <-   colnames(obj@data@data)
             colnames(obj@param@mu) <-  paste("class-",1:length(obj@param@pi),sep="")
             colnames(obj@param@sd) <-  paste("class-",1:length(obj@param@pi),sep="")
             proba <- rep(0, obj@data@n)
             for (j in 1:obj@data@d){
               obj@param@mu[j,1] <- mean(obj@data@data[which(obj@data@notNA[,j]==1),j])
               obj@param@sd[j,1] <- sd(obj@data@data[which(obj@data@notNA[,j]==1),j])
               proba[which(obj@data@notNA[,j]==1)] <- proba[which(obj@data@notNA[,j]==1)] + dnorm(obj@data@data[which(obj@data@notNA[,j]==1),j], obj@param@mu[j,1], obj@param@sd[j,1], log = TRUE)
             }
             obj@partitions@zMAP <- rep(1, obj@data@n)
             obj@partitions@zOPT <- rep(1, obj@data@n)
             obj@partitions@tik <- matrix(1, obj@data@n, 1)
             obj@criteria@loglikelihood <- sum(proba)
             obj@criteria@BIC <- obj@criteria@loglikelihood - obj@data@d*log(obj@data@n)
             obj@criteria@ICL <- ICLcontinuous(obj) 
             obj@criteria@MICL <- obj@criteria@ICL
             obj@criteria@degeneracyrate <- 0
             return(obj)         
           }
)