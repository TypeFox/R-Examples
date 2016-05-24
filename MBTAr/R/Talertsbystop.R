Talertsbystop <- function(stop_id,include_access_alerts=FALSE,include_service_alerts=TRUE,api_key){
  query <- "alertsbystop"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&stop=",stop_id,"&include_access_alerts=",include_access_alerts,"&include_service_alerts=",include_service_alerts,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- NULL
  for(i in 1:length(dl$alerts$alert_id)){
    alltimes <- NULL
    for(j in 1:nrow(dl$alerts$effect_periods[[i]])){
      thistime <- data.frame(
        stop_id = dl$stop_id,
        stop_name = dl$stop_name,
        alert_id = dl$alerts$alert_id[i],
        effect_name = dl$alerts$effect_name[i],
        effect = dl$alerts$effect[i],
        cause= dl$alerts$cause[i],
        header_text = dl$alerts$header_text[i],
        short_header_text = dl$alerts$short_header_text[i],
        description_text = dl$alerts$description_text[i],
        severity = dl$alerts$severity[i],
        created_dt = dl$alerts$created_dt[i],
        last_modified_dt = dl$alerts$last_modified_dt[i],
        service_effect_text = dl$alerts$service_effect_text[i],
        timeframe_text = dl$alerts$timeframe_text[i],
        alert_lifecycle = dl$alerts$alert_lifecycle[i],
        effect_start = dl$alerts$effect_periods[[i]]$effect_start[j],
        effect_end = dl$alerts$effect_periods[[i]]$effect_end[j],
        affected_route_type = ifelse(length(dl$alerts$affected_services$services[[i]]$route_type)>0,dl$alerts$affected_services$services[[i]]$route_type,NA),
        affected_mode_name = ifelse(length(dl$alerts$affected_services$services[[i]]$mode_name)>0,dl$alerts$affected_services$services[[i]]$mode_name,NA),
        affected_route_id = ifelse(length(dl$alerts$affected_services$services[[i]]$route_id)>0,dl$alerts$affected_services$services[[i]]$route_id,NA),
        affected_route_name = ifelse(length(dl$alerts$affected_services$services[[i]]$route_name)>0,dl$alerts$affected_services$services[[i]]$route_name,NA),
        affected_direction_id = ifelse(length(dl$alerts$affected_services$services[[i]]$direction_id)>0,dl$alerts$affected_services$services[[i]]$direction_id,NA), # sometimes NA
        affected_direction_name = ifelse(length(dl$alerts$affected_services$services[[i]]$direction_name)>0,dl$alerts$affected_services$services[[i]]$direction_name,NA), # sometimes NA
        affected_trip_id = ifelse(length(dl$alerts$affected_services$services[[i]]$trip_id)>0,dl$alerts$affected_services$services[[i]]$trip_id,NA), # sometimes NA
        affected_trip_name = ifelse(length(dl$alerts$affected_services$services[[i]]$trip_name)>0,dl$alerts$affected_services$services[[i]]$trip_name,NA), # sometimes NA
        affected_stop_id = ifelse(length(dl$alerts$affected_services$services[[i]]$stop_id)>0,dl$alerts$affected_services$services[[i]]$stop_id,NA), # sometimes NA
        affected_stop_name = ifelse(length(dl$alerts$affected_services$services[[i]]$stop_name)>0,dl$alerts$affected_services$services[[i]]$stop_name,NA), # sometimes NA
        affected_route_hide = ifelse(length(dl$alerts$affected_services$services[[i]]$route_hide)>0,dl$alerts$affected_services$services[[i]]$route_hide,NA), # sometimes NA
        affected_elev_id = ifelse(length(dl$alerts$affected_services$elevators[[i]]$elev_id)>0,dl$alerts$affected_services$elevators[[i]]$elev_id,NA),
        affected_elev_name = ifelse(length(dl$alerts$affected_services$elevators[[i]]$elev_name)>0,dl$alerts$affected_services$elevators[[i]]$elev_name,NA),
        affected_elev_type = ifelse(length(dl$alerts$affected_services$elevators[[i]]$elev_type)>0,dl$alerts$affected_services$elevators[[i]]$elev_type,NA),
        affected_elev_stop_id = ifelse(length(dl$alerts$affected_services$elevators[[i]]$stop$stop_id)>0,dl$alerts$affected_services$elevators[[i]]$stop$stop_id,NA),
        affected_elev_stop_name = ifelse(length(dl$alerts$affected_services$elevators[[i]]$stop$stop_name)>0,dl$alerts$affected_services$elevators[[i]]$stop$stop_name,NA),
        affected_elev_stop_parent_id = ifelse(length(dl$alerts$affected_services$elevators[[i]]$stop$parent_station)>0,dl$alerts$affected_services$elevators[[i]]$stop$parent_station,NA),
        affected_elev_stop_parent_name = ifelse(length(dl$alerts$affected_services$elevators[[i]]$stop$parent_station_name)>0,dl$alerts$affected_services$elevators[[i]]$stop$parent_station_name,NA)
      )
      alltimes <- rbind(alltimes,thistime)
    }
    allout <- rbind(allout,alltimes)
  }
  return(allout)
}
