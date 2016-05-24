hclust.wrap <-
function(x,centers) {
   hc <- hclust(dist(x))

   res=cutree(hc, k = centers) #k = 1 is trivial
   return(res)
}
