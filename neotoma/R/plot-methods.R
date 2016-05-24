#' @importFrom graphics legend par plot points
#' @export
plot.site <- function(x, database = NULL, ...) {
  if(is.null(database)) {
    plot(x$long, x$lat, xlab = "Longitude", ylab = "Latitude", ...)
  }  else {
    plot(x$long, x$lat, xlab = "Longitude", ylab = "Latitude", type = 'n')
    points(x$long, x$lat, ...)
    
  }
}

#' @export
plot.dataset <- function(x, ...) {
    site <- get_site(x)
    plot(site, ...)
}

#' @export
plot.dataset_list <- function(x, ...) {
  site <- get_site(x)
  
  types <- factor(sapply(x, function(x)x$dataset.meta$dataset.type))
  
  levels(types) <- gsub(' ', '\n', levels(types))
  
  old_par <- par()
  par(mar=c(5, 4, 4, 5))
  plot(site, pch=as.numeric(types)-1, ..., bty='L')
  legend(max(site$long)+0.5,max(site$lat), levels(types), 
         pch=(1:length(types))-1, xpd=TRUE, cex = 0.75)
  par(mar = old_par$mar)
}

#' @export
plot.download <- function(x, ...) {
  site <- get_site(x)
  plot(site, ...)
}

#' @export
plot.download_list <- function(x, ...) {
  site <- get_site(x)

  types <- factor(sapply(x, function(x)x$dataset$dataset.meta$dataset.type))
  levels(types) <- gsub(' ', '\n', levels(types))
  
  old_par <- par()
  par(mar=c(5, 4, 4, 5))
  plot(site, pch=as.numeric(types)-1, ..., bty='L')
  legend(max(site$long)+0.5,max(site$lat), levels(types), 
         pch=(1:length(types))-1, xpd=TRUE, cex = 0.75)
  par(mar = old_par$mar)

}
