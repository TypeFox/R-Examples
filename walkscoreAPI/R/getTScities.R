getTScities <-
function(key){
     URL <- paste("http://transit.walkscore.com/transit/supported/cities/?wsapikey=",key,sep="")
     X <- character(0)
     X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))

     citylines <- grep("\"city\":",X)
     statelines <- citylines + 1

     citylist <- X[citylines]
     citylist <- gsub("\"city\":","",citylist)
     citylist <- gsub("  ","",citylist)
     citylist <- gsub("\"","",citylist,fixed=TRUE)
     citylist <- gsub(",","",citylist,fixed=TRUE)

     statelist <- X[statelines]
     statelist <- gsub("\"state\":","",statelist)
     statelist <- gsub("  ","",statelist)
     statelist <- gsub("\"","",statelist,fixed=TRUE)
     statelist <- gsub(",","",statelist,fixed=TRUE)
 
     print(paste(citylist,statelist))
}

