tv_episode_images <-
function(api_key, id, season_number, episode_number){
    
    url <- fromJSON(GET(url=paste("http://api.themoviedb.org/3/tv/", id, "/season/", season_number, 
                                  "/episode/", episode_number, "/images?api_key=", 
                                  api_key, sep=""))$url)
    
    return(url)
    
}
