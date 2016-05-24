#' @export
print.download <- function(x, ...) {

  site <- x$dataset$site.data$site.name
  
  if (all(is.na(x$sample.meta[,c('age.older', 'age', 'age.younger')]))) {
    age.set <- c(NA, NA)
  } else {
    age.set <- suppressWarnings(range(as.vector(x$sample.meta[,c('age.older', 'age', 'age.younger')]), 
                                      na.rm = TRUE))
  }
  
  age.set[!is.finite(age.set)] <- NA
  names(age.set) <- c('age.younger', 'age.older')
  
  locs <- as.numeric(get_site(x)[,c('long', 'lat')])
  
  types <- get_dataset(x)[[1]]$dataset.meta$dataset.type
  if (!is.na(x$dataset$access.date)) {
    cat(paste0('A download object for ',
             x$dataset$site.data$site.name, '\n',
             'Accessed ', format(x$dataset$access.date, "%Y-%m-%d %H:%M"), 'h. \n'))
  } else {
    cat(paste0('A download object for ',
               x$dataset$site.data$site.name, '\n',
               'No associated access date. \n'))
  }
  print(format(data.frame(dataset.id = x$dataset$dataset.meta$dataset.id, 
                          site.name = site, 
                          long = locs[1],
                          lat = locs[2],
                          age.young = age.set[1],
                          age.old = age.set[2],
                          type = types), 
               justify='left'), row.names=FALSE)
  
  NULL
}
