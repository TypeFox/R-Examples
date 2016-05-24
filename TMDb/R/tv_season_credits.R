tv_season_credits <-
function(api_key, id, season_number){
    
    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/season/", season_number, 
                                  "/credits?api_key=", api_key, sep=""))$url)
    
    return(url)
    
}
