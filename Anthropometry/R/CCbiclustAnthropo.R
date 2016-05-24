CCbiclustAnthropo <- function(data,waistVariable,waistCirc,lowerVar,nsizes,nBic,diffRanges,percDisac,dir){
  
 n <- c()
 dims <- list() ; res <- list() ; wr <- list() ; mat <- list() ; tab_acc <- list() ; 
 ColBics <- list() ; delta <- list() ; disac <- list()
  
 setwd(dir)
  
 for (i in 1 : (nsizes-1)){ 
  da_size <- data[(waistVariable >= waistCirc[i]) & (waistVariable < waistCirc[i + 1]), ] 
  n[i] <- dim(da_size)[1]
    
  da_size2 <- da_size[,lowerVar] 
  da_size3 <- da_size2 / 10 
    
  diff_ranges <- as.vector(apply(da_size3, 2, range)[2,] - apply(da_size3, 2, range)[1,])
  selectedCols <- which(diff_ranges >= diffRanges[[i]][1] & diff_ranges <= diffRanges[[i]][2])
      
  da_size4 <- as.matrix(da_size3[,selectedCols])
  dims[[i]] <- dim(da_size4)
      
  delta[[i]] <- 1
  disac[[i]] <- nrow(da_size4)
      
  while (disac[[i]] > ceiling(percDisac * nrow(da_size4))){
        
   res[[i]] <- biclust(da_size4, method = BCCC(), delta = delta[[i]], alpha = 1.5, number = nBic[i]) 
   wr[[i]] <- writeclust(res[[i]], row = TRUE, noC = res[[i]]@Number)
   disac[[i]] <- sum(wr[[i]] == 0) 
   
   if(disac[[i]] > ceiling(percDisac * nrow(da_size4))){
    delta[[i]] = delta[[i]] + 1
   }
  }
    
  if(res[[i]]@Number == 0){
   mat[[i]] <- NA
  }else{
    mat[[i]] <- sapply(1 : res[[i]]@Number, overlapBiclustersByRows, res[[i]])
    acc <- c()
    for(j in 1:nrow(mat[[i]])){
     acc[j] <- length(which(mat[[i]][j,] == 0)) 
    }
      
    tab_acc[[i]] <- list(table(acc), nBic[i], dims[[i]][1], disac[[i]])
    rm(acc)
     
    tab_acc[[i]] <- list(table(acc), nBic[i], dims[[i]][1])
    rm(acc)
      
    fils <- res[[i]]@RowxNumber
    cols <- res[[i]]@NumberxCol
   
    fun1 <- function(x){
     list(rownames(da_size4[fils[, x], ]), colnames(da_size4[, cols[x,]])) 
    }
      
    nom_var <- paste("list.bicl", i, sep = ".")
    assign(nom_var, lapply(1:length(fils[1, ]), fun1)) 
      
    lista <- get(paste("list.bicl.", i, sep = ""))
    save(lista, file = paste("list.bicl", i, "RData", sep = "."))
      
    filsBic <- list()
    colsBic <- list()
    for(j in 1:res[[i]]@Number){
     filsBic[[j]] <- biclusternumber(res[[i]], number = 1:res[[i]]@Number)[[j]][1] 
     colsBic[[j]] <- biclusternumber(res[[i]], number = 1:res[[i]]@Number)[[j]][2]
    }
    tab <- table(unlist(colsBic)) 
    dimnames_tab <- attr(tab, "dimnames")
    dimnames_tab2 <- as.numeric(unlist(dimnames_tab))
    ColBics[[i]] <- colnames(da_size2)[dimnames_tab2[tab == nBic[i]]] 
   }
 }
  
 return(list(res=res,dims=dims,delta=delta,disac=disac,mat=mat,tab_acc=tab_acc,ColBics=ColBics))
}