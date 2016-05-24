#' @export
print.dataset_list <- function(x, ...){
  
  if(length(x) == 0){
    cat('Empty dataset list.')
    return()
  }
  
  dates <- range(sapply(x, '[[', 'access.date'))
  sites <- sapply(lapply(x, '[[', 'site.data'), '[[', 'site.name')
  dataset.id <- sapply(lapply(x, '[[', 'dataset.meta'), '[[', 'dataset.id')
  types <- sapply(lapply(x, '[[', 'dataset.meta'), '[[', 'dataset.type')
  
  #  Get site locations:
  locs <- get_site(x)[,c('long', 'lat')]
  
  cat(paste0('A dataset_list containing ', length(x), ' objects:\n',
           'Accessed from ', 
           format(as.POSIXct(dates[1], origin=Sys.time()-as.numeric(Sys.time())), "%Y-%m-%d %H:%M"),
           'h to ',
           format(as.POSIXct(dates[2], origin=Sys.time()-as.numeric(Sys.time())), "%Y-%m-%d %H:%M"),
           'h. \n',
           'Datasets:\n'))
  print(format(data.frame(dataset.id, 
                          site.name = sites, 
                          long= locs[,1],
                          lat = locs[,2],
                          type = types), 
               justify='left'), row.names=FALSE)
  
  NULL
}
