`medv` <- 
function(psm, h=0.99){
    cutree(hclust(as.dist(1-psm),method="complete"), h=h)
}
