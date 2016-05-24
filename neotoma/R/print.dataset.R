#' @export
print.dataset <- function(x, ...){

  site <- as.character(x$site$site.name)
  
  locs <- as.numeric(get_site(x)[,c('long', 'lat')])
  types <- x$dataset.meta$dataset.type
  
  if(!is.na(x$access.date)){
    cat(paste0('A dataset for ',
             x$site$site.name, '\n',
             'Accessed ', format(x$access.date, "%Y-%m-%d %H:%M"), 'h. \n'))
  } else {
    cat(paste0('A dataset for ',
               x$site$site.name), '\n')
  }
  
  print(format(data.frame(dataset.id = x$dataset.meta$dataset.id, 
                          site.name = site, 
                          long = locs[1],
                          lat = locs[2],
                          type = types),
               justify='left'), row.names=FALSE)
  
  NULL
}
