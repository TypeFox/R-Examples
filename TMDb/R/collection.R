collection <-
function(api_key, id, language=NA, append_to_response=NA){
    
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(append_to_response) && !(append_to_response %in% "images")){
        stop("append_to_response can be NA or images")
    }
        
    l <- list(language=language, append_to_response=append_to_response)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/collection/", id, "?api_key=", 
                                      api_key, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/collection/", id, "?api_key=", 
                                      api_key, sep=""))$url)        
    }
    
    return(url)
   
}
