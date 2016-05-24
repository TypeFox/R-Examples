company <-
function(api_key, id, append_to_response=NA){
    
    if(!is.na(append_to_response) && !(append_to_response %in% "movies")){
        stop("append_to_response can be NA or movies string")
    }
    
    l <- list(append_to_response=append_to_response)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/company/", id, "?api_key=", 
                                      api_key, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/company/", id, "?api_key=", 
                                      api_key, sep=""))$url)        
    }
    
    return(url)
    
}
