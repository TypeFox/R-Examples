plotRule <-
function(Tree) {

  if (length(Tree[["id"]]) > 1) {
    goWhile <- TRUE
  
    while(goWhile){
  
      # max density of each cluster 
      Tops <- numeric(length(Tree[["id"]]))
      for (j in Tree[["id"]]) {
        Tops[j] <- max(Tree[["density"]][Tree[["DataPoints"]][[j]]])
      }
    
      
      uniqueParents <- unique(Tree$parent)
      uniqueParentsNo0 <- setdiff(uniqueParents, 0)
      
      
      for (i in uniqueParentsNo0) {
        
        bros <- Tree[["sons"]][[i]]
        newID <- bros
        orderID <- bros[order(Tops[bros], decreasing = TRUE)]
        newID[order(Tops[bros], decreasing=TRUE)] <- bros  

        # I have corrected the way NewID are assigned
        # TODO what if same height?
            
        if (sum(bros != newID) != 0) {
        
          NewTree <- Tree
          
          #update new IDs
          NewTree[["lambdaTop"]][bros] <- Tree[["lambdaTop"]][orderID]  
          NewTree[["rTop"]][bros] <- Tree[["rTop"]][orderID]  
          NewTree[["kappaTop"]][bros] <- Tree[["kappaTop"]][orderID]  
          NewTree[["alphaTop"]][bros] <- Tree[["alphaTop"]][orderID]  
          
          for (j in seq(along = bros)) {
            NewTree[["parent"]][which(Tree[["parent"]] == bros[j])] <- newID[j]
          }
            
          for (j in seq(along = newID)) {
            if (!is.null(Tree[["sons"]][bros[j]][[1]])){
              NewTree[["sons"]][[newID[j]]] <- Tree[["sons"]][bros[j]][[1]]
            } else {
              NewTree[["sons"]][[newID[j]]] <- NA
            }
          }
          
          for (j in seq(along = newID)) {
            NewTree[["DataPoints"]][[newID[j]]] <- Tree[["DataPoints"]][[bros[j]]]
          }
          
          ## Now we modify Xbase and silos
          for (s in seq(along = NewTree[["id"]])) {
            if (NewTree$parent[s] == 0) {
            Bros <- which(NewTree[["parent"]] == 0) 
            rank <- which(Bros == s)
            NewTree[["silo"]][[s]] <- siloF(c(0, 1), length(Bros), rank)  
            NewTree[["Xbase"]][s] <- sum(NewTree[["silo"]][[s]]) / 2
            } else {
            Bros <- which(NewTree[["parent"]] == NewTree[["parent"]][s])  
            rank <- which(Bros == s)
            NewTree[["silo"]][[s]] <- siloF(NewTree[["silo"]]
                [[NewTree[["parent"]][s]]], length(Bros), rank)
            NewTree[["Xbase"]][s] <- sum(NewTree[["silo"]][[s]]) / 2      
            }
          }
          
          Tree <- NewTree 
          break
        }     
      }
      
      if (i == rev(uniqueParentsNo0)[1]) {
        goWhile <- FALSE
      }    
    }
  }
  return(Tree)
}
