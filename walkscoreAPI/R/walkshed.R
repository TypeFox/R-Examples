walkshed <- function(x,y,key){
     URL <- paste("http://api.walkscore.com/walk_shed?lat=",y,"&lon=",x,"&wsapikey=",key,sep="")
     X <- character(0)
     X <- c(X, scan(file = URL, what = "", sep = " ", quiet = TRUE))

     status <- X[grep("{status:",X,fixed=TRUE) + 1]
     status <- as.numeric(gsub(",","",status,fixed=TRUE))
     lat <- X[grep("{lat",X,fixed=TRUE) + 1]
     lat <- as.numeric(gsub(",","",lat))
     lon <- X[grep("{lat",X,fixed=TRUE) + 3]
     lon <- as.numeric(gsub("},","",lon))
     geo <- X[grep("{type:",X,fixed=TRUE) + 1]
     geo <- gsub(",","",geo)
     cstart <- grep("coordinates:",X,fixed=TRUE) + 1
     cend <- grep("radius:",X,fixed=TRUE) - 1
     lonlist <- c()
     latlist <- c()
     count <- 0
     for (i in cstart:cend){
         count <- count + 1
             if (count %% 2 == 1){
                n <- X[i]
                n <- gsub("[","",n,fixed=TRUE)
                n <- gsub(",","",n,fixed=TRUE)
                n <- as.numeric(n)
                lonlist <- c(lonlist,n)
             }
             if (count %% 2 == 0){
                n <- X[i]
                n <- gsub("]","",n,fixed=TRUE)
                n <- gsub(",","",n,fixed=TRUE)
                n <- gsub("}","",n,fixed=TRUE)
                n <- as.numeric(n)
                latlist <- c(latlist,n)
             }
         }
     coords <- data.frame(lonlist,latlist)
     rad <- X[grep("radius:",X,fixed=TRUE) + 1]
     rad <- as.numeric(gsub("},","",rad,fixed=TRUE))
     slat <- X[grep("snapped_lat:",X,fixed=TRUE) + 1]
     slat <- as.numeric(gsub(",","",slat,fixed=TRUE))
     slon <- X[grep("snapped_lon:",X,fixed=TRUE) + 1]
     slon <- as.numeric(gsub("}","",slon,fixed=TRUE))

     object <- list()
     class(object) <- "Walkshed"
     object$status <- status
     object$origin <- c(lon,lat)
     object$geometry <- geo
     object$coordinates <- coords
     object$radius <- rad 
     object$snappedLong <- slon
     object$snappedLat <- slat
     return(object)
}