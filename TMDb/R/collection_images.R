collection_images <-
function(api_key, id, language=NA, append_to_response=NA, include_image_language=NA){
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(append_to_response) && !(append_to_response %in% "images")){
        stop("append_to_response can be NA or images string")
    }
    
    if((!is.na(include_image_language) && !is.character(include_image_language)) 
       || (length(include_image_language)>5)){
        stop("include_image_language must be comma separated string with max 5 elements ISO 69-1 code.")
    }
    
    l <- list(language=language, append_to_response=append_to_response, 
              include_image_language=include_image_language)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/collection/", id, "/images?api_key=", 
                                      api_key, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/collection/", id, "/images?api_key=", 
                                      api_key, sep=""))$url)        
    }
    
    return(url)
    
}
