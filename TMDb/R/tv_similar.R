tv_similar <-
function(api_key, id, page=1, language=NA, append_to_response=NA){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(append_to_response)){
        append_to_response<-gsub(pattern = " ",replacement = "", x = append_to_response)
        v<-unlist(strsplit(append_to_response, split=","))
        tv_method <- c("alternative_titles", "changes", "content_ratings", "credits", "external_ids",
        "images", "keywords", "similar", "translations", "videos", "latest", "on_the_air",
        "airing_today", "top_rated", "popular")
        for (i in v){
            if (!(i %in% tv_method))
                stop(paste(i,  "is not a valid tv_method"))
        }
    }
    
    l <- list(append_to_response=append_to_response)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/similar?api_key=", 
                                      api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/similar?api_key=", 
                                      api_key, sep=""))$url)
    }
    
    return(url)
    
}
