dendro.subjects <-
function(data){  
#function(data, type=list()){  
# !! to be done: allow also asymmetric binary variables
  
  D.subjects <- dist.subjects(data)
  #D.subjects <- dist.subjects(data, type=type)
  as.dendrogram(hclust(D.subjects))
}
