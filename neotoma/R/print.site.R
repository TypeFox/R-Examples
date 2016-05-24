#' @export

print.site <- function(x, ...){
  class(x) <- 'data.frame'
  
  print(format(data.frame(site.name = x$site.name, 
                          long = x$long,
                          lat = x$lat,
                          elev = x$elev),
               justify='left'), row.names=FALSE)

  cat(paste0('A site object containing ',nrow(x),' sites and 8 parameters.\n'))
  
}
