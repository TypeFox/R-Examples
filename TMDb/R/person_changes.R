person_changes <-
function(api_key, id, start_date=NA, end_date=NA){
    
    l <- list(start_date=start_date, end_date=end_date)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/", id, "/changes?api_key=", 
                                      api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/", id, "/changes?api_key=", 
                                      api_key, sep=""))$url)
    }
    
    return(url)
    
}
