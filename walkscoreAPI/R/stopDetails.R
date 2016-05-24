stopDetails <-
function(stopid,key){
     URL <- paste("http://transit.walkscore.com/transit/stop/",stopid,"/?wsapikey=",key,sep="")
     X <- character(0)
     X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
 
     string <- X[grep("\"name\":",X,fixed=TRUE)]
     string2 <- X[grep("\"lon\":",X,fixed=TRUE)]
     string3 <- X[grep("\"lat\":",X,fixed=TRUE)]
     rbegin <- grep("\"route_ids\":",X,fixed=TRUE)
     rend <- grep("],",X,fixed=TRUE)
 
    if (length(X) > 0){
       name <- strsplit(string,": ")
       name <- gsub("\"","",name[[1]][2])
       lon <- strsplit(string2,": ")
       lon <- gsub(",","",lon[[1]][2])
       lon <- as.numeric(lon)
       lat <- strsplit(string3,": ")
       lat <- gsub(",","",lat[[1]][2])
       lat <- as.numeric(lat)
       
       rlist <- c()
       for (i in (rbegin+1):(rend-1)){
            str <- X[i]
            str <- gsub(" ","",str)
            str <- gsub("\"","",str,fixed=TRUE)
            str <- gsub(",","",str)
            rlist <- c(rlist,str)
          }
       }

    else {
       print("Error, invalid stop ID")
    }
    object <- list()
    class(object) <- "StopDetails"
    object$stopID <- stopid
    object$stopName <- name
    object$stopLong <- lon
    object$stopLat <- lat
    object$routelist <- rlist
    return(object)
  }

