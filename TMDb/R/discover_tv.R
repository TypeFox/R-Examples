discover_tv <-
function(api_key, page=1, language=NA, sort_by=NA, first_air_date_year=NA,
                        vote_count.gte=NA, vote_average.gte=NA, with_genres=NA,
                        with_networks=NA, first_air_date.gte=NA, first_air_date.lte=NA){
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.na(first_air_date_year) && !is.character(first_air_date_year)){
        stop("first_air_date_year must be a date string like YYYY")
    }
    
    if(!is.na(first_air_date.gte) && !is.character(first_air_date.gte)){
        stop("first_air_date.gte must be a date string like YYYY-MM-DD")
    }
    
    if(!is.na(first_air_date.lte) && !is.character(first_air_date.lte)){
        stop("first_air_date.lte must be a date string like YYYY-MM-DD")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(!is.na(vote_count.gte) && !is.wholenumber(vote_count.gte)){
        stop("vote_count.gte must be an integer.")
    }
    
    if(!is.na(vote_average.gte) && !is.wholenumber(vote_average.gte)){
        stop("vote_average.gte must be a float.")
    }
    
    l <- list(language=language, sort_by=sort_by, first_air_date_year=first_air_date_year,
              vote_count.gte=vote_count.gte, vote_average.gte=vote_average.gte, 
              with_genres=with_genres, with_networks=with_networks, first_air_date.gte=first_air_date.gte, 
              first_air_date.lte=first_air_date.lte)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/discover/tv?api_key=", 
                                      api_key, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/discover/tv?api_key=", 
                                      api_key, sep=""))$url)        
    }
    
    return(url)
    
}
