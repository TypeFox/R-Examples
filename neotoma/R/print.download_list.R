#' @export
print.download_list <- function(x, ...){
    
  if(any(sapply(x, is.null))){
    if(all(sapply(x, is.null))){
      stop('All returned datasets are empty.\n')
    } else{
      x <- x[!sapply(x, is.null)]
      empty <- sum(sapply(x, is.null))
    }
  } else {
    empty <- 0
  }
  
  
  dates <- range(sapply(lapply(x, '[[', 'dataset'), '[[', 'access.date'))
  sites <- sapply(lapply(lapply(x, '[[', 'dataset'), '[[', 'site.data'), '[[', 'site.name')
  dataset.id <- sapply(lapply(lapply(x, '[[', 'dataset'), '[[', 'dataset.meta'), '[[', 'dataset.id')
  
  types <- sapply(lapply(get_dataset(x), '[[', 'dataset.meta'), '[[', 'dataset.type')
  
  #  Get minimum and maximum ages from the dataset object within a download:
  site_ages <- function(x){
    if(all(is.na(x$sample.meta[,c('age.older', 'age', 'age.younger')]))){
      age.set <- c(NA, NA)
    } else{
      age.set <- suppressWarnings(range(as.vector(x$sample.meta[,c('age.older', 'age', 'age.younger')]), na.rm=TRUE))
    }
    
    age.set
    
  }
  
  age.set <- suppressWarnings(t(sapply(x, site_ages)))
  
  age.set[!is.finite(age.set)] <- NA
  colnames(age.set) <- c('age.younger', 'age.older')
  
  #  Get site locations:
  locs <- get_site(x)[,c('long', 'lat')]

  cat(paste0('A download_list containing ', length(x), ' objects:\n',
           'Accessed from ', 
           format(as.POSIXct(dates[1], origin=Sys.time()-as.numeric(Sys.time())), "%Y-%m-%d %H:%M"),
           'h to ',
           format(as.POSIXct(dates[2], origin=Sys.time()-as.numeric(Sys.time())), "%Y-%m-%d %H:%M"),
           'h. \n',
           'Datasets:\n'))
  print(format(data.frame(dataset.id, 
                          site.name = sites, 
                          locs, 
                          age.set,
                          type = types), 
               justify='left'), row.names=FALSE)
  
  if(empty == 1){ cat('There is one empty download associated with this download_list\n')}
  if(empty > 1){ cat(paste0('There are ',empty, ' downloads associated with this download_list\n'))}
  
  NULL
}
