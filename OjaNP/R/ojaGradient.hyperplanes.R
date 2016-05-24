`ojaGradient.hyperplanes` <-
function(Hyperplanes,x){
   B <- apply(Hyperplanes,1,function(y){ojaGradient.hyperplane(d=y,x=x)})
   if (is.matrix(B))
      return(rowMeans(B))
   else
      return(mean(B))
}
