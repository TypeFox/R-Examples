Talertbyid <- function(alert_id,include_access_alerts=FALSE,include_service_alerts=TRUE,api_key){
  query <- "alertbyid"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&id=",alert_id,"&include_access_alerts=",include_access_alerts,"&include_service_alerts=",include_service_alerts,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- NULL
  for(i in 1:length(dl$alert_id)){
    alltimes <- NULL
    for(j in 1:nrow(dl$effect_periods)){
      thistime <- data.frame(
        alert_id = dl$alert_id[i],
        effect_name = dl$effect_name[i],
        effect = dl$effect[i],
        cause= dl$cause[i],
        header_text = dl$header_text[i],
        short_header_text = dl$short_header_text[i],
        description_text = dl$description_text[i],
        severity = dl$severity[i],
        created_dt = dl$created_dt[i],
        last_modified_dt = dl$last_modified_dt[i],
        service_effect_text = dl$service_effect_text[i],
        timeframe_text = dl$timeframe_text[i],
        alert_lifecycle = dl$alert_lifecycle[i],
        effect_start = dl$effect_periods$effect_start[j],
        effect_end = dl$effect_periods$effect_end[j],
        affected_route_type = ifelse(length(dl$affected_services$services$route_type)>0,dl$affected_services$services$route_type,NA),
        affected_mode_name = ifelse(length(dl$affected_services$services$mode_name)>0,dl$affected_services$services$mode_name,NA),
        affected_route_id = ifelse(length(dl$affected_services$services$route_id)>0,dl$affected_services$services$route_id,NA),
        affected_route_name = ifelse(length(dl$affected_services$services$route_name)>0,dl$affected_services$services$route_name,NA),
        affected_direction_id = ifelse(length(dl$affected_services$services$direction_id)>0,dl$affected_services$services$direction_id,NA), # sometimes NA
        affected_direction_name = ifelse(length(dl$affected_services$services$direction_name)>0,dl$affected_services$services$direction_name,NA), # sometimes NA
        affected_trip_id = ifelse(length(dl$affected_services$services$trip_id)>0,dl$affected_services$services$trip_id,NA), # sometimes NA
        affected_trip_name = ifelse(length(dl$affected_services$services$trip_name)>0,dl$affected_services$services$trip_name,NA), # sometimes NA
        affected_stop_id = ifelse(length(dl$affected_services$services$stop_id)>0,dl$affected_services$services$stop_id,NA), # sometimes NA
        affected_stop_name = ifelse(length(dl$affected_services$services$stop_name)>0,dl$affected_services$services$stop_name,NA), # sometimes NA
        affected_route_hide = ifelse(length(dl$affected_services$services$route_hide)>0,dl$affected_services$services$route_hide,NA), # sometimes NA
        affected_elev_id = ifelse(length(dl$affected_services$elevators$elev_id)>0,dl$affected_services$elevators$elev_id,NA),
        affected_elev_name = ifelse(length(dl$affected_services$elevators$elev_name)>0,dl$affected_services$elevators$elev_name,NA),
        affected_elev_type = ifelse(length(dl$affected_services$elevators$elev_type)>0,dl$affected_services$elevators$elev_type,NA),
        affected_elev_stop_id = ifelse(length(dl$affected_services$elevators$stop$stop_id)>0,dl$affected_services$elevators$stop$stop_id,NA),
        affected_elev_stop_name = ifelse(length(dl$affected_services$elevators$stop$stop_name)>0,dl$affected_services$elevators$stop$stop_name,NA),
        affected_elev_stop_parent_id = ifelse(length(dl$affected_services$elevators$stop$parent_station)>0,dl$affected_services$elevators$stop$parent_station,NA),
        affected_elev_stop_parent_name = ifelse(length(dl$affected_services$elevators$stop$parent_station_name)>0,dl$affected_services$elevators$stop$parent_station_name,NA)
      )
      alltimes <- rbind(alltimes,thistime)
    }
    allout <- rbind(allout,alltimes)
  }
  return(allout)
}
