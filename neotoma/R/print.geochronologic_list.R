#' @export
print.geochronologic_list <- function(x, ...){
  
  dates <- range(sapply(x, function(y)y[[1]]$access.date))
  sites <- sapply(lapply(x, function(y)y[[1]]$site.data), '[[', 'site.name')
  dataset.id <- sapply(lapply(x, function(y)y[[1]]$dataset.meta), '[[', 'dataset.id')
  
  #  Get site locations:
  locs <- do.call(rbind.data.frame,
                  lapply(x, function(y)get_site(y[[1]])[,c('long', 'lat')]))
  
  ages <- sapply(x, function(y) nrow(y[[2]]))
  mins <- sapply(x, function(y) min(y[[2]]$age, na.rm = TRUE))
  maxs <- sapply(x, function(y) max(y[[2]]$age, na.rm = TRUE))
  
  cat(paste0('A geochronology_list containing ', length(x), ' objects:\n',
           'Accessed from ', 
           format(as.POSIXct(dates[1], origin=Sys.time()-as.numeric(Sys.time())), "%Y-%m-%d %H:%M"),
           'h to ',
           format(as.POSIXct(dates[2], origin=Sys.time()-as.numeric(Sys.time())), "%Y-%m-%d %H:%M"),
           'h. \n',
           'Geochronologies:\n'))
  print(format(data.frame(id = dataset.id, 
                          site.name = sites, 
                          long= locs[,1],
                          lat = locs[,2],
                          ages = ages,
                          min  = mins,
                          max  = maxs,
                          interval = (maxs - mins) / ages),
               justify='left'), row.names=FALSE)
  
  NULL
}
