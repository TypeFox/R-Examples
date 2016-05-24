########################################################################################################################
## La fonction DesignOutput permet de mettre en forme les parametres en fonction de leur nature VSLCMresultsContinuous
## ou VSLCMresultsCategorical. Elle est appelee a la fin de l'estimation des parametres
########################################################################################################################
setGeneric ( name= "DesignOutput",  def = function(reference){ standardGeneric("DesignOutput")})

## Cas de variables continues
setMethod( f = "DesignOutput", 
           signature(reference="VSLCMresultsContinuous"), 
           definition = function(reference){
             reference@model@omega <-  as.numeric(reference@model@omega)
             names(reference@model@omega) <- colnames(reference@data@data)
             if (reference@model@g>1){
               if (reference@strategy@vbleSelec==FALSE){
                 reference@partitions@zOPT <- numeric()
                 reference@criteria@MICL <- numeric()
                 reference@partitions@zMAP <- as.numeric(reference@partitions@zMAP) + 1
               }
               # On remet les valeurs manquantes
               for (j in 1:reference@data@d){
                 if (any(reference@data@notNA[,j]==0))
                   reference@data@data[which(reference@data@notNA[,j]==0),j] <- NA
               }
               if (reference@strategy@paramEstim){
                 if (reference@criteria@degeneracyrate != 1){
                   rownames(reference@param@mu)  <-   colnames(reference@data@data)
                   rownames(reference@param@sd)  <-   colnames(reference@data@data)
                   if (reference@criteria@degeneracyrate>0.5)
                     warning(paste("The rate of degeneracy for the EM algorithm is", reference@criteria@degeneracyrate ), call. = FALSE)
                   reference@param@pi <- as.numeric(reference@param@pi)
                   names(reference@param@pi) <-  paste("class-",1:length(reference@param@pi),sep="")
                   colnames(reference@param@mu) <-  paste("class-",1:length(reference@param@pi),sep="")
                   colnames(reference@param@sd) <-  paste("class-",1:length(reference@param@pi),sep="")
                   
                   reference@partitions@zMAP <- as.numeric(reference@partitions@zMAP) + 1
                   reference@partitions@zOPT <- as.numeric(reference@partitions@zOPT) + 1
                   colnames(reference@partitions@tik) <-  paste("class-",1:reference@model@g,sep="")
                   reference@criteria@BIC <- reference@criteria@loglikelihood  - 0.5*(reference@model@g-1 + reference@model@g*2*sum(reference@model@omega) + 2*sum(1-reference@model@omega))*log(reference@data@n)
                   reference@criteria@ICL <- ICLcontinuous(reference) 
                 }else{
                   warning("All the models got error (degeneracy)", call. = FALSE)
                 }
               }
             }
             return(reference)
           }
)