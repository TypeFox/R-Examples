person_movie_credits <-
function(api_key, id, language=NA, append_to_response=NA){
    
    if(!is.na(append_to_response)){
        append_to_response<-gsub(pattern = " ",replacement = "", x = append_to_response)
        v<-unlist(strsplit(append_to_response, split=","))
        people_method <- c("movie_credits", "tv_credits", "combined_credits",
        "external_ids", "images", "tagged_images", "changes", "popular", "latest")
        for (i in v){
            if (!(i %in% people_method))
                stop(paste(i,  "is not a valid people_method"))
        }
    }
    
    l <- list(language=language, append_to_response=append_to_response)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/", id, "/movie_credits?api_key=", 
                                      api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/person/", id, "/movie_credits?api_key=", 
                                      api_key, sep=""))$url)
    }
    
    return(url)
    
}
