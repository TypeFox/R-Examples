check_lake <- function(lake_name){
  
  if (!tolower(lake_name) %in% tolower(get_lake_names())){
    stop(paste0('lake_name=', lake_name, ' not recognized in GLTC dataset'))
  }
}

check_climate <- function(climate_name){
  in_db <- tolower(climate_name) %in% tolower(get_climate_names())
  if (any(!in_db)){
    stop(paste0('climate_name=', climate_name[!in_db], ' not recognized in GLTC dataset'))
  }
}