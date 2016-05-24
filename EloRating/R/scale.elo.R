#scale.elo <- function(x) 
scale.elo <- function (x, center = TRUE, scale = TRUE) {round((x-min(x, na.rm=T))/max((x-min(x, na.rm=T)), na.rm=T), 3)}
