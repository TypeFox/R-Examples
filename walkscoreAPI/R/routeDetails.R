routeDetails <-
function(routeid,key){
     URL <- paste("http://transit.walkscore.com/transit/route/",routeid,"/?wsapikey=",key,sep="")
     X <- character(0)
     X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
 
     string <- X[grep("\"name\":",X,fixed=TRUE)]
     string2 <- X[grep("\"category\":",X,fixed=TRUE)]
     string3 <- X[grep("\"agency\":",X,fixed=TRUE)]
     string4 <- X[grep("\"agency_url\":",X,fixed=TRUE)]
     string5 <- X[grep("\"geometry_wkt\":",X,fixed=TRUE)]
     sbegin <- grep("\"stop_ids\":",X,fixed=TRUE)
     send <- grep("],",X,fixed=TRUE)
 
    if (length(X) > 0){
       name <- strsplit(string,": ")
       name <- gsub("\"","",name[[1]][2])
       name <- gsub(", ","",name)
       cat <- strsplit(string2,": ")
       cat <- gsub("\"","",cat[[1]][2], fixed=TRUE)
       cat <- gsub(", ","",cat)
       age <- strsplit(string3,": ")
       age <- gsub("\"","",age[[1]][2],fixed=TRUE)
       age <- gsub(", ","",age)
       url <- strsplit(string4,": ")
       url <- gsub("\"","",url[[1]][2],fixed=TRUE)
       url <- gsub(", ","",url)
       geom <- strsplit(string5,": ")
       geom <- gsub("\"","",geom[[1]][2],fixed=TRUE)
       geom <- gsub(", ","",geom)
       geom <- gsub("LINESTRING","",geom)
       geom <- gsub("(","",geom,fixed=TRUE)
       geom <- gsub(")","",geom,fixed=TRUE)
       geom <- strsplit(geom,",")
       if (length(geom)>1){
           geolist <- data.frame()
           for (i in 1:length(geom[[1]])){
               coords <- strsplit(geom[[1]][i]," ")
               geolist[i,1] <- as.numeric(coords[[1]][1])
               geolist[i,2] <- as.numeric(coords[[1]][2])
           }
       names(geolist) <- c("X","Y")
       }
       else {
           geolist <- "Unavailable" 
       }      
       slist <- c()
       for (i in (sbegin+1):(send-1)){
            str <- X[i]
            str <- gsub(" ","",str)
            str <- gsub("\"","",str,fixed=TRUE)
            str <- gsub(",","",str)
            slist <- c(slist,str)
          }
       }

    else {
       print("Error, invalid route ID")
    }
    object <- list()
    class(object) <- "RouteDetails"
    object$routeID <- routeid
    object$routeName <- name
    object$routeCatagory <- cat
    object$agengy <- age
    object$agencyURL <- url 
    object$routeGeometry <- geolist
    object$stopList <- slist
    return(object)
  }

