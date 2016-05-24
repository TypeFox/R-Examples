#' Removes missing values (rows and column to obtain a large full matrix)
#' @export
#' @description Removes missing values (rows and column alternatively) to obtain a large full matrix
#' @param X the dataset (matrix) with missing values
#' @param nbstep number of cutting steps (may remove several rows or columns at each step)
#' @param std (boolean) remove constant covariates
#' @param verbose (boolean) to print the result
#' @param coercing vector of the covariates to keep (names or index)
#' @param Xout (boolean) to export or not the reduced matrix (if not, indices are sufficient)
#' @return 
#' \item{individus_restants}{Index of remaining individuals}
#' \item{variables_restantes}{Index of remaining variables}
#' \item{X}{If Xout=TRUE, the reduced dataset without missing values}

#' @examples
#'    \dontrun{

#'    data<-mtcars
#'    require(CorReg)
#'   datamiss=Terminator(target = data,wrath=0.05)#5% of missing values
#'   datamiss
#'   showdata(datamiss)#plot positions of the missing values
#'   reduced=Conan(X=datamiss)
#'   reduced
#'     }


Conan<-function(X=X,nbstep=Inf,std=FALSE,verbose=FALSE,coercing=NULL,Xout=TRUE){
   steploc=1 
   nloc=nrow(X)
   ploc=ncol(X)
   n=nloc
   p=ploc
   individus_restants=1:n
   variables_restantes=1:p
   if(!is.null(coercing)){#on a des variables a garder  
      #on va purger les individus manquants sur ces variables (pas le choix)
      for(i in 1:length(coercing)){
         if(is.numeric(coercing[i])){
            j=coercing[i]
         }else{
            j=which(names(X)==coercing[i])[1] #si doublon on prend le premier
         }
         if(anyDuplicated(c(1:ploc,j))>0){#si on ne tape pas a cote
            quimank=which(is.na(X[,j]))
            if(length(quimank)>0){
               X=X[-quimank,]
               individus_restants=individus_restants[-quimank]
            }
         }
      }
      nloc=nrow(X)
      ploc=ncol(X)
   }
   #nbmank/prod(dim(basefull))
   if(nbstep==Inf){
      nbstep=max(c(nloc,ploc))
   }
   if(verbose){print(dim(X))}
   if(std){#nettoyage des constantes si demande
      sd_loc<-function(vect){
         return(sd(vect[!is.na(vect)]))
      }
      quiconst=which(apply(X,2,sd_loc)==0)
      if(length(quiconst)>0){
         X=X[,-quiconst]
         variables_restantes=variables_restantes[-quiconst]      
      }   
   }
   if(verbose){print(dim(X))}
   quimank=which(is.na(X),arr.ind = TRUE)
   nbmank=dim(quimank)[1]
   M=Matrix::Matrix(0,nrow = nrow(X),ncol=ncol(X))
   M[quimank]=1 
   while (steploc<=nbstep & nbmank>0 & nloc>0 & ploc>0){
      if(verbose){print(steploc)}     
      nbmank_var=colSums(M)/nloc#taux=nbmank colonne/ nombre de ligne
      nbmank_ind=rowSums(M)/ploc
      if(max(nbmank_var)>max(nbmank_ind)){#on supprime variable
         quisuppr=which(nbmank_var==max(nbmank_var))
         X=X[,-quisuppr]
         M=M[,-quisuppr]
         variables_restantes=variables_restantes[-quisuppr]
      }else{#on supprime ligne
         if(verbose){print("lign")}
         quisuppr=which(nbmank_ind==max(nbmank_ind))
         if(length(quisuppr)>0){#si il y a quelque chose a supprimer
            X=X[-quisuppr,]
            M=M[-quisuppr,]
            individus_restants=individus_restants[-quisuppr]
            if(verbose){print(dim(X))}
            if(std){#nettoyage des constantes si demande
               quiconst=which(apply(X,2,sd_loc)==0)
               if(length(quiconst)>0){
                  X=X[,-quiconst]
                  variables_restantes=variables_restantes[-quiconst]
               }        
            }
         }
      }
      if(verbose){print(dim(X)) }
      nloc=nrow(X)
      ploc=ncol(X)
      quimank=which(is.na(X),arr.ind = TRUE)
      nbmank=dim(quimank)[1]
      steploc=steploc+1
   }
   if(verbose){
      print(paste("I've killed ",n-nloc,"individuals and ",p-ploc," covariates."))
   }
   if(Xout){
      return(list(X=X,individus_restants=individus_restants,variables_restantes=variables_restantes))
   }else{
      return(list(individus_restants=individus_restants,variables_restantes=variables_restantes))
   }
}