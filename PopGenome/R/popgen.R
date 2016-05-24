popgen <- function(Code_matrix,Populations=FALSE,outgroup=FALSE,methods=FALSE,include.unknown=TRUE,gff=FALSE,FAST,SNP.DATA){


## KONVERTING ##########################################
# Check if an allel doesnt exists in the alignement ####
# ----------------------------------------------### ####
# If it doesnt exist delete it #########################

ALLROWS <- FALSE

if(!is.list(Populations)){                     #### wenn keine Population definiert nimm alle !
   Populations  <- list(1:dim(Code_matrix)[1])
   Populations2 <- list(rownames(Code_matrix))
   npops        <- 1
   ALLROWS      <- TRUE    

}else{

  npops        <- length(Populations)
  Populations2 <- vector("list",npops) # beinhaltet die Namen der Sequenzen

for(xx in 1:npops){
  
   if(is.character(Populations[[xx]])){
  
     namesX  <- Populations[[xx]]
     isit    <- which(is.na(match(namesX,rownames(Code_matrix))))
   
     if(length(isit)>=1){
   # Sequence doesnt exist
        if(length(attr(Code_matrix,"path"))!=0){
           warning("GEN",attr(Code_matrix,"path"))
        }
      warning("The allele ----->: ",namesX[isit], " <---- doesnt exist")
      namesX <- namesX[-isit]
     }
   
   # unbekannte geloescht ! 
    Populations[[xx]]   <- match(namesX,rownames(Code_matrix))
    Populations2[[xx]]  <- namesX
  } # End Population is character
 
  if(is.numeric(Populations[[xx]])){
     Populations2[[xx]] <- rownames(Code_matrix)[Populations[[xx]]]
  }
  
 }# End of for npops
}# else

# --------------------------------------------------------------------------------#

########################################################
# Check if Population == NULL # Verkleiner sonst die Population
########################################################

if(is.list(Populations)){
   res          <- delNULLpop(Populations)
   Populations  <- res$Populations
   popmissing   <- res$popmissing
   res          <- delNULLpop(Populations2)
   Populations2 <- res$Populations
   npops        <- npops - length(popmissing)
}
#########################################################



#---------------------------------------------------------------------------------#
# Outgroup #################### --------------------------------------------------#

outgroup2 <- vector()        # outgroup Namen

if(is.character(outgroup)){

 namesX  <- outgroup
 isit    <- which(is.na(match(namesX,rownames(Code_matrix))))
 
 if(length(isit)>=1){

   if(length(attr(Code_matrix,"path"))!=0){
      warning("region.stats",attr(Code_matrix,"path"))
   }
  warning("The outgroup allele ----->: ",namesX[isit], "<---- doesnt exist")
  namesX   <- namesX[-isit] # unbekannte geloescht !
 
 } 
  outgroup  <- match(namesX,rownames(Code_matrix))
  outgroup2 <- rownames(Code_matrix)[outgroup]
  
} # end of outgroup character

if(is.numeric(outgroup)){
  outgroup2 <- rownames(Code_matrix)[outgroup]
}


#####################################
#### NAMES ##########################
#####################################

NAME    <- paste("pop",1:npops)
if(outgroup[1]){
 NAME <- c(NAME,"outgroup")
}


#dat   <- new("DATA") # create a new class of type DATA
   
     # Save genename # ----------------------------------
     if(length(attr(Code_matrix,"path"))==0){genename <- "unknown"}else{
      genename <- attr(Code_matrix,"path")}
     # ---------------------------------------------------

   # GET INFORMATION FROM THE CODING MATRIX    <---------------------------------------------- take only the defined population alleles
   if(outgroup[1]==FALSE){outgr <- NULL}else{outgr <- outgroup}
   #
   if(!ALLROWS){Code_matrix <- Code_matrix[unique(c(unlist(Populations),outgr)),,drop=FALSE]}
   #
   if(dim(Code_matrix)[1]<=1){return(NA)} 
   # --------------------------------------
   
   # Aktualisiere Populationen auf geschrumpfte Code_matrix
  
   if(!ALLROWS){          # Wenn nicht alle Sequencen (beziehungweise keine Population definiert wurde)
     for(xx in 1:npops){
        Populations[[xx]] <- match(Populations2[[xx]],rownames(Code_matrix))
     } 
   }

   if(outgroup[1]!=FALSE){outgroup <- match(outgroup2,rownames(Code_matrix))}
   if(length(outgroup)==0){outgroup <- FALSE} # Outgroup existiert garnicht !
   ###### GET DATA #############################


   obj <- get_data(Code_matrix,include.unknown,gff=gff,FAST,SNP.DATA)
   

   ## Exception
   if(length(obj)==1){ # No statistics calculated
   return(NA) # no biallelic sites
   }
   ############

   if(length(rownames(Code_matrix)) > 0){ 
    genes       <- rownames(Code_matrix)
   }else{
    genes       <- 1:dim(Code_matrix)[1]
   }

   #dat@MATRIX      <- Code_matrix
   
   #if(is.list(Populations)){
   populations  <- Populations
   populations2 <- Populations2
   popmissing   <- popmissing
   #}
  
   outgroup                       <- as.matrix(outgroup)   # bezieht sich auf das komplette Alignment
   outgroup2                      <- as.matrix(outgroup2)
   rownames(obj$biallelic.matrix) <- genes
   
   #if(!is.list(Populations)){return(dat)}


  return(list(data.sum=obj,genename=genename,genes=genes,populations=populations,populations2=populations2,popmissing=popmissing,outgroup=outgroup,outgroup2=outgroup2))

}
