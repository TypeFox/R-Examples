stopSearch <-
function(x,y,key){
   URL <- paste("http://transit.walkscore.com/transit/search/stops/?lat=",y,"&lon=",x,"&wsapikey=",key,sep="")  
   X <- character(0)
   X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
   split <- strsplit(X,"{",fixed=TRUE)
   stoplines <- grep("distance",split[[1]])
   stoplines <- c(stoplines,length(split[[1]])+1)

   if (length(X) > 0){
       slist <- list()
       for (i in 1:(length(stoplines)-1)){
           text <- split[[1]][stoplines[i]:(stoplines[(i+1)]-1)]
           info <- text[grep("distance",text)] 
           info <- strsplit(info, ",", fixed = TRUE)
           dist <- info[[1]][grep("distance",info[[1]])]
           dist <- strsplit(dist,": ")
           dist <- as.numeric(dist[[1]][2])
           sname <- info[[1]][grep("name",info[[1]])]
           sname <- strsplit(sname,": ",fixed=TRUE)
           sname <- gsub("\"","",sname[[1]][2],fixed=TRUE)
           lastline <- text[length(text)]
           lastline <- strsplit(lastline,", ", fixed = TRUE)
           lon <- lastline[[1]][grep("\"lon",lastline[[1]],fixed=TRUE)]
           lon <- strsplit(lon,": ",fixed=TRUE)
           lon <- as.numeric(lon[[1]][2])
           lat <- lastline[[1]][grep("\"lat",lastline[[1]],fixed=TRUE)]
           lat <- strsplit(lat,": ",fixed=TRUE)
           lat <- as.numeric(lat[[1]][2])
           sid <- lastline[[1]][length(lastline[[1]])]
           sid <- strsplit(sid,": ",fixed=TRUE)
           sid <- gsub("\"","",sid[[1]][2],fixed=TRUE)
           sid <- gsub("}","",sid,fixed=TRUE)
           sid <- gsub("]","",sid,fixed=TRUE)

           rlist <- list()
           r <- grep("category",text)

           for (j in 1:length(r)){
               rtext <- strsplit(text[r[j]],",")
               cat <- rtext[[1]][grep("category",rtext[[1]])]
               cat <- strsplit(cat,": ",fixed=TRUE)
               cat <- gsub("\"","",cat[[1]][2],fixed=TRUE)
               id <- rtext[[1]][grep("\"id\"",rtext[[1]],fixed=TRUE)]
               id <- strsplit(id,": ",fixed=TRUE)
               id <- gsub("\"","",id[[1]][2],fixed=TRUE)
               rname <- rtext[[1]][grep("\"name\"",rtext[[1]],fixed=TRUE)]
               rname <- strsplit(rname,": ",fixed=TRUE)
               rname <- gsub("\"","",rname[[1]][2],fixed=TRUE)
               rname <- gsub("}","",rname,fixed=TRUE)
               rname <- gsub("]","",rname,fixed=TRUE)
               age <- rtext[[1]][grep("\"agency\"",rtext[[1]],fixed=TRUE)]
               age <- strsplit(age,": ",fixed=TRUE)
               age <- gsub("\"","",age[[1]][2],fixed=TRUE)
               url <- rtext[[1]][grep("\"agency_url\"",rtext[[1]],fixed=TRUE)]
               url <- strsplit(url,": ",fixed=TRUE)
               url <- gsub("\"","",url[[1]][2],fixed=TRUE)
               rlist <- list()
               class(rlist) <- "RouteDetails"
               rlist$id <- id
               rlist$routeName <- rname
               rlist$routeCatagory <- cat
               rlist$routeAgency <- age
               rlist$routeURL <- url
          }
        object <- list()
        class(object) <- "Stop2"
        object$stopName <- sname
        object$stopID <- id
        object$stopDistance <- dist
        object$stopLong <- lon
        object$stopLat <- lat
        object$routeDetails <- rlist
        slist[[i]] <- object
       }
   }
 return(slist)
}

