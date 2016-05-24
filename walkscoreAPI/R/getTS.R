getTS <-
function(x,y,city,state,key){
    city <- gsub(" ", "+",city)
    URL <- paste("http://transit.walkscore.com/transit/score/?lat=",y,"&lon=",x,"&city=",city,"&state=",state,"&wsapikey=",key,sep="")
    X <- character(0)
    X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))

    string <- X[grep("transit_score",X)]
    string2 <- X[grep("ws_link",X)]
    string3 <- X[grep("description",X)]
    string4 <- X[grep("summary",X)]

    if (length(X) > 0){
       tscore <- strsplit(string,": ")
       tscore <- gsub(",","",tscore[[1]][2])
       tscore <- as.numeric(tscore)
       link <- strsplit(string2,": ")
       link <- gsub("\"","",link[[1]][2],fixed=TRUE)
       link <- gsub(", ","",link)
       desc <- strsplit(string3,": ")
       desc <- gsub("\"","",desc[[1]][2],fixed=TRUE)
       desc <- gsub(", ","",desc)
       sum <- strsplit(string4,": ")
       sum <- gsub("\"","",sum[[1]][2],fixed=TRUE)
       sum <- gsub(", ","",sum)
       }
   else {
       print("error, please check supported cities list")
  }

  object <- list()
  class(object) <- "TransitScore"
  object$transitscore <- tscore
  object$url <- link
  object$description <- desc
  object$summary <- sum 
  return(object)
}

