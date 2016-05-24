#cutree.graph
"cuttree.graph"<-function(tree){
if (inherits(tree,"hclust")){
  plot(tree, hang=-1, main="Choice of the number of clusters by cutting the dendrogram")
  coupe<-locator(n=1)
  while (coupe$y<0){ 
    cat("no cluster defined \n")
    coupe<-locator(n=1)
  }
  #abline(h=coupe$y, col=2)
  rect.hclust(tree, h=coupe$y)
} 
else stop("the tree must be 'hclust' or 'agnes' class")       
classes<-cutree(as.hclust(tree),h=coupe$y)
return(classes)
}
