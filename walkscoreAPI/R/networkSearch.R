networkSearch <-
function(x,y,key){
   URL <- paste("http://transit.walkscore.com/transit/search/network/?lat=",y,"&lon=",x,"&wsapikey=",key,sep="")  
   X <- character(0)
   X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
   split <- strsplit(X,"}",fixed=TRUE)
   routesline <- grep("\"routes\"",split[[1]],fixed=TRUE)
   stopsline <- grep("\"stops\"",split[[1]],fixed=TRUE)
   rtext <- split[[1]][routesline:(stopsline-1)]
   stext <- split[[1]][stopsline:length(split[[1]])]
   rtext <- gsub("{\"routes\":","",rtext,fixed=TRUE)
   stext <- gsub("\"stops\":","",stext,fixed=TRUE)   
   routelist <- list()
   stoplist <- list()

   if (length(X) > 0){
       for (i in 1:(length(rtext)-1)){
           text <- strsplit(rtext[i],", ",fixed=TRUE)
           cat <- text[[1]][grep("category",text[[1]])]
           cat <- strsplit(cat,": ",fixed=TRUE)
           cat <- gsub("\"","",cat[[1]][3],fixed=TRUE)
           age <- text[[1]][grep("\"agency\"",text[[1]],fixed=TRUE)]
           age <- strsplit(age,": ",fixed=TRUE)
           age <- gsub("\"","",age[[1]][2],fixed=TRUE)
           name <- text[[1]][grep("\"name\"",text[[1]],fixed=TRUE)]
           name <- strsplit(name,": ",fixed=TRUE)
           name <- gsub("\"","",name[[1]][2],fixed=TRUE)
           url <- text[[1]][grep("\"agency_url\"",text[[1]],fixed=TRUE)]
           url <- strsplit(url,": ",fixed=TRUE)
           url <- gsub("\"","",url[[1]][2],fixed=TRUE)
           id <- text[[1]][grep("\"id\"",text[[1]],fixed=TRUE)]
           id <- strsplit(id,": ",fixed=TRUE)
           id <- gsub("\"","",id[[1]][2],fixed=TRUE)
           sstart <- grep("\"stop_ids\"",text[[1]])
           send <- grep("\"id\"",text[[1]])-1
           slist <- c()
               for (j in sstart:send){
                   stop <- gsub("\"stop_ids\": [","",text[[1]][j],fixed=TRUE)
                   stop <- gsub("\"","",stop,fixed=TRUE)
                   stop <- gsub("]","",stop,fixed=TRUE)
                   slist <- c(slist,stop)
               }
           object <- list()
           class(object) <- "Route"
           object$routeID <- id
           object$routeName <- name
           object$routeCatagory <- cat
           object$agency <- age
           object$agencyURL <- url
           object$stopList <- slist
           routelist[[i]] <- object
          }
       
       for (i in 1:(length(stext)-2)){
           text <- strsplit(stext[i],", ",fixed=TRUE)
           lat <- text[[1]][grep("\"lat\"",text[[1]],fixed=TRUE)]
           lat <- strsplit(lat,": ",fixed=TRUE)
           lat <- as.numeric(lat[[1]][3])
           lon <- text[[1]][grep("\"lon\"",text[[1]],fixed=TRUE)]
           lon <- strsplit(lon,": ",fixed=TRUE)
           lon <- as.numeric(lon[[1]][2])
           name <- text[[1]][grep("\"name\"",text[[1]],fixed=TRUE)]
           name <- strsplit(name,": ",fixed=TRUE)
           name <- gsub("\"","",name[[1]][2],fixed=TRUE)
           id <- text[[1]][grep("\"id\"",text[[1]],fixed=TRUE)]
           id <- strsplit(id,": ",fixed=TRUE)
           id <- gsub("\"","",id[[1]][2],fixed=TRUE)
           rstart <- grep("\"route_ids\"",text[[1]])
           rend <- grep("\"lon\"",text[[1]])-1
           rlist <- c()
               for (j in rstart:rend){
                   route <- gsub("\"route_ids\": [","",text[[1]][j],fixed=TRUE)
                   route <- gsub("\"","",route,fixed=TRUE)
                   route <- gsub("]","",route,fixed=TRUE)
                   rlist <- c(rlist,route)
               }
           object <- list()
           class(object) <- "Stop"
           object$stopID <- id
           object$stopName <- name
           object$stopLong <- lon
           object$stopLat <- lat
           object$routeList <- rlist
           stoplist[[i]] <- object
          }
      }
     
      object <- list()
      class(object) <- "NetworkSearch"
      object$routelist <- routelist
      object$stoplist <- stoplist
      return(object)      
}

