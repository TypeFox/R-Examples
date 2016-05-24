#' @export
print.geochronologic <- function(x, ...){

  site <- as.character(x[[1]]$site$site.name)
  
  locs <- as.numeric(get_site(x[[1]])[,c('long', 'lat')])
  
  cat(paste0('Geochronological data for ',
           x[[1]]$site$site.name, '\n',
           'Accessed ', format(x[[1]]$access.date, "%Y-%m-%d %H:%M"), 'h. \n'))

  if(nrow(x[[2]])>1){
    interval <- mean(diff(x[[2]]$age))
  } else {
    interval <- NA
  }
  
  print(format(data.frame(dataset.id = x[[1]]$dataset.meta$dataset.id, 
                          site.name = site, 
                          long = locs[1],
                          lat = locs[2],
                          ages = nrow(x[[2]]),
                          min  = min(x[[2]]$age, na.rm = TRUE),
                          max  = max(x[[2]]$age, na.rm = TRUE),
                          interval = interval),
               justify='left'), row.names=FALSE)
  
  NULL
}
