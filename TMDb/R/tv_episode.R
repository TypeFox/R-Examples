tv_episode <-
function(api_key, id, season_number, episode_number, language=NA, append_to_response=NA){
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(append_to_response)){
        append_to_response<-gsub(pattern = " ",replacement = "", x = append_to_response)
        v<-unlist(strsplit(append_to_response, split=","))
        tv_episode_method <- c("changes", "credits", "external_ids", "images", "videos")
        for (i in v){
            if (!(i %in% tv_episode_method))
                stop(paste(i,  "is not a valid tv_episode_method"))
        }
    }
    
    l <- list(language=language, append_to_response=append_to_response)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/season/", season_number, "/episode/",
                                      episode_number, "?api_key=", api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/season/", season_number, "/episode/",
                                      episode_number, "?api_key=", api_key, sep=""))$url)
    }
    
    return(url)
    
}
