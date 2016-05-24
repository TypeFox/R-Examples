sort_silhouette <-
function(sil,cluster)  { 
     ss <- NULL 
     ks <- sort(unique(cluster)) 
     for(k in ks) {
          ss <- c(ss,sort(sil[cluster==k]))
     }
     ss
}
