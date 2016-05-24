#'
#' @title
#' Indices values and Jack-knife indices for a single topology.
#'
#' @description
#' The funtion calculates standard and terminal jack-knifed indices I and W 
#' [see Miranda-Esquivel 2016], along with Posadas et al. 2001 modifications.
#' 
#' @param tree is a single tree with n terminals, an ape phylo object.
#' 
#' @param distribution species distributions in n areas, a data.frame
#' 
#' @param jtip is the proportion of terminals to delete, real (range 0-1).
#' 
#' @param verbose Boolean. If TRUE, the output reports the number of deleted terminals/topologies. 
#' 
#' @param standard "distribution" or "tree" to standarize by the 
#' by the sum of indices in the distribution or  the sum of indices in the tree
#' 
#' @examples
#' library(jrich)
#' data(tree)
#' data(distribution)
#' ## Standarized by the sum of indices in the distribution
#' Calculate.Index(tree=tree, distribution = distribution, verbose=TRUE, standard = "distribution")
#' 
#' ## Standarized by the sum of indices in the tree (as figure 1 in Miranda-Esquivel 2016)
#' Calculate.Index(tree=tree, distribution = distribution, verbose=TRUE, standard = "tree")
#' 
#'
#'@author Miranda-Esquivel Daniel R.
#'
#'





Calculate.Index <- function (tree = tree, distribution = distribution, jtip = 0, verbose = TRUE, standard = "distribution") {

  ## Errors on trees / distributions
  ## names and numbers
  
    if (names(distribution)[1]=="especie"){
    names(distribution)[1] <- "species"
    }
  
    if (length(tree$tip.label) == length(distribution$species)){
        if (all(tree$tip.label[order(tree$tip.label)] == 
                distribution$species[order(distribution$species)])){
      }else{
          stop("distributions and tree(s) MUST have the same names for species and terminals")
         } 
  }else{
      stop("distributions and tree(s) MUST have the same number of species and terminals")    
       }
   

  	areas       <-  names(distribution)[-length(distribution)]
  	
  	 if(length(areas) < 2) {warning("Endemism values could be missleading")}
  	
  	especies    <-  distribution$species
  
  	deleted.Terminals    <-  0
  
  	for (i in 1:length(especies)){
    
    	if (jtip > runif(1)) {
        	#! print(paste("especie ",especies[i]," sera borrada",sep=""))
      	deleted.Terminals    <-  deleted.Terminals + 1
            
        distribution[i,-(length(areas)+1)] <- rep (0,length(areas))
    	}    
  	} 
  
  	if(verbose){
    print(paste("Deleted",deleted.Terminals,"out of",length(especies),sep="  "))
  	}
  
  	filas<-length(names(distribution))-1
  	resultados <-as.data.frame(matrix(data=0,nrow=filas,ncol=13))
  	names(resultados)<-c("area","I","Ie","Is","Ise","W","We","Ws","Wse","rich","endem","jtopol","jtip")
  	resultados$area <- names(distribution)[names(distribution)!="species"]

  
  ##
  ## si arboles es multiphylo hacer el calculo por arbol, sumar los indices y promediar
  ## como los arboles no tienen (o si?) la misma secuencia de terminales
  ## ordenarlos por terminales
  ##

  	tree <- reorder.phylo(tree,order="postorder")
  
  	W <- IndexW(tree=tree)  
  	I <- IndexI(tree=tree)
  
  	names(I) <- names(W) <- tree$tip.label
  
  	match <- match(tree$tip.label,distribution$species)
  
  	distribution <- distribution[match,]
  	distribution <- distribution[,names(distribution)!="species"]

## Omar Leon reported
## 2016 - 08 - 25
## At this point, the numeric variables are treated as characters,
## so the next operations can't be made, then I treat all numeric
## variables as numeric, finally the calculations work.

        for (i in 1:length(colnames(distribution))){
           distribution[,i] <- as.numeric(distribution[,i])
        }
  	
  	## here I have to include a possible solution to handle a single area 
  
    resultados$rich <-  apply(as.matrix(distribution),2,sum)
  
  
  endemicSpecies       <-   apply(as.matrix(distribution),1,sum)
  endemicSpecies[which(endemicSpecies != 1)] = 0
  endemicityMatrix     <-  endemicSpecies*distribution
  sumEndemicityMatrix  <-  apply(as.matrix(endemicityMatrix),2,sum) 
  resultados$endem     <-  resultados$rich*sumEndemicityMatrix
  
  #if (resultados$rich < resultados$endem ){
  #    resultados$endem=resultados$rich
  #    }
  
  	indiceI.areas <- I*distribution
  	indiceW.areas <- W*distribution
  

  	resultados$I <-  apply(as.matrix(indiceI.areas),2,sum)
  	resultados$W <-  apply(as.matrix(indiceW.areas),2,sum)

  
    resultados$Ie <-  resultados$I/resultados$rich
    resultados$We <-  resultados$W/resultados$rich

  	##
  	## Two different approaches to s, per area or per Index,
  	## note that in both cases the proportions are the same while
  	## the absolute values differ
  	##
  	if (standard == "distribution"){
      
  		resultados$Is <- resultados$I/sum(resultados$I)
  		resultados$Ws <- resultados$W/sum(resultados$W)	
      
  		resultados$Ise <- resultados$Ie/sum(resultados$Ie)  
  		resultados$Wse <- resultados$We/sum(resultados$We)  
      
  		}else{
  		
        resultados$Is <- resultados$I/sum(I)
  			resultados$Ws <- resultados$W/sum(W)
        
  			resultados$Ise <- resultados$Ie/sum(I)
  			resultados$Wse <- resultados$We/sum(W)
  		}
  

  

 
  
  	resultados[,c(2:13)] <- round(resultados[,c(2:13)],3)
  
  #! if(jtip<0){jtip<-0}
  resultados$jtip <- jtip
  
  resultados[resultados=="NaN"] <- 0
  
  return(resultados)
}
