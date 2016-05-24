silhouette.km <-
function(x,centers) { 
     dd <- NULL
     k <- nrow(centers)
     for(i in 1:k) {
          xr <- sweep(x,2,centers[i,],'-')
          dd<-cbind(dd,apply(xr^2,1,sum))
     }
     dd <- apply(dd,1,sort)[1:2,]
     (dd[2,]-dd[1,])/dd[2,]
}
