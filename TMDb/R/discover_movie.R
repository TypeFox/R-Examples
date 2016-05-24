discover_movie <-
function(api_key, certification_country=NA, certification=NA, certification.lte=NA, 
                           include_adult=FALSE, include_video=TRUE, language=NA, page=1,
                           primary_release_year=NA, primary_release_date.gte=NA, primary_release_date.lte=NA,
                           release_date.gte=NA, release_date.lte=NA, sort_by=NA, vote_count.gte=NA,
                           vote_count.lte=NA, vote_average.gte=NA, vote_average.lte=NA, with_cast=NA,
                           with_crew=NA, with_companies=NA, with_genres=NA, with_keywords=NA,
                           with_people=NA, year=NA){

    if(!is.na(certification_country) && !is.character(certification_country)){
        stop("certification_country must be a ISO3166-1 code")
    }
    
    if(!is.na(certification_country) && is.na(certification.lte)){
        stop("certification.lte is required when certification_country is specified.")
    }
    
    if(!is.logical(include_adult)){
        stop("include_adult must have a logical value.")
    }
    
    if(!is.logical(include_video)){
        stop("include_video must have a logical value.")
    }
    
    if(!is.na(language) && !is.character(language)){
        stop("language must be a ISO639-1 code")
    }
    
    if(page<1 || page>1000){
        stop("page must be a number between 1 and 1000")
    }
    
    if(!is.na(primary_release_year) && !is.character(primary_release_year)){
        stop("primary_release_year must be a date string like YYYY")
    }
    
    if(!is.na(primary_release_date.gte) && !is.character(primary_release_date.gte)){
        stop("primary_release_date.gte must be a date string like YYYY-MM-DD")
    }
    
    if(!is.na(primary_release_date.lte) && !is.character(primary_release_date.lte)){
        stop("primary_release_date.lte must be a date string like YYYY-MM-DD")
    }
    
    if(!is.na(release_date.gte) && !is.character(release_date.gte)){
        stop("release_date.gte must be a date string like YYYY-MM-DD")
    }
    
    if(!is.na(release_date.lte) && !is.character(release_date.lte)){
        stop("release_date.lte must be a date string like YYYY-MM-DD")
    }
    
    if(!is.na(vote_count.gte) && !is.wholenumber(vote_count.gte)){
        stop("vote_count.gte must be an integer.")
    }
    
    if(!is.na(vote_count.lte) && !is.wholenumber(vote_count.lte)){
        stop("vote_count.lte must be an integer.")
    }
    
    if(!is.na(vote_average.gte) && !is.wholenumber(vote_average.gte)){
        stop("vote_average.gte must be a float.")
    }
    
    if(!is.na(vote_average.lte) && !is.wholenumber(vote_average.lte)){
        stop("vote_average.lte must be a float.")
    }
    
    if(!is.na(year) && !is.character(year)){
        stop("year must be a date string like YYYY")
    }
    
    l <- list(certification_country=certification_country, certification=certification, 
              certification.lte=certification.lte, include_adult=include_adult, 
              include_video=include_video, language=language, page=page, 
              primary_release_year=primary_release_year, primary_release_date.gte=primary_release_date.gte, 
              primary_release_date.lte=primary_release_date.lte, release_date.gte=release_date.gte, 
              release_date.lte=release_date.lte, sort_by=sort_by, vote_count.gte=vote_count.gte,
              vote_count.lte=vote_count.lte, vote_average.gte=vote_average.gte, 
              vote_average.lte=vote_average.lte, with_cast=with_cast, with_crew=with_crew, 
              with_companies=with_companies, with_genres=with_genres, with_keywords=with_keywords,
              with_people=with_people, year=year)
    l <- l[!is.na(l)]
    
    if(length(l)>0){
        params <- paste("&", names(l), "=",l, sep="", collapse="")
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/discover/movie?api_key=", 
                                      api_key, params, sep=""))$url)   
    } else{
        url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/discover/movie?api_key=", 
                                      api_key, sep=""))$url)        
    }
    
    return(url)
    
}
