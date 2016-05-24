checkTSsupport <-
function(city,state,key){
     ci <- tolower(city)
     st <- tolower(state)
     URL <- paste("http://transit.walkscore.com/transit/supported/cities/?wsapikey=",key,sep="")
     X <- character(0)
     X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
     
     ciline <- 0
     ciline <- c(grep(ci,X),ciline)
   
     if (ciline[1] > 0){
        stline <- 0 
        stline <- c(grep(st,X[ciline[1]+1]),stline)
             if (stline[1] > 0){
                 return(TRUE)
             }
             else {
                return(FALSE)
             }
       }
     else {
        return(FALSE)
     }
}

