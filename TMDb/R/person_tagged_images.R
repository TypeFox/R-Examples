person_tagged_images <-
function(api_key, id, page=1, language=NA){
    
    l <- list(language=language, page=page)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/", id, "/tagged_images?api_key=", 
                                      api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/", id, "/tagged_images?api_key=", 
                                      api_key, sep=""))$url)
    }
    
    return(url)
    
}
