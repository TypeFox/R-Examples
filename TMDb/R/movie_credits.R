movie_credits <-
function(api_key, id, append_to_response=NA){
    
    if(!is.na(append_to_response)){
        append_to_response<-gsub(pattern = " ",replacement = "", x = append_to_response)
        v<-unlist(strsplit(append_to_response, split=","))
        movie_method <- c("alternative_titles", "credits", "images", "keywords", "releases", "videos", "translations", "similar", "reviews", "lists", "changes", "rating", "latest", "upcoming", "now_playing", "popular", "top_rated")
        for (i in v){
            if (!(i %in% movie_method))
                stop(paste(i,  "is not a valid movie_method"))
        }
    }
    
    l <- list(append_to_response=append_to_response)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/movie/", id, "/credits?api_key=", 
                                      api_key, params, sep=""))$url)
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/movie/", id, "/credits?api_key=", 
                                      api_key, sep=""))$url)
    }
    
    return(url)
    
}
