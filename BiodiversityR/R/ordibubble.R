`ordibubble` <-
function(ordiplot,var,...) {
    y2 <- abs(var)
    ordiscores <- scores(ordiplot, display="sites")
    for (i in 1:length(var)) {
        if (var[i] < 0) {
            var[i] <- NA
        }else{
            y2[i] <- NA
        }
    }
    if (sum(var,na.rm=T) > 0) {graphics::symbols(ordiscores[,1], ordiscores[,2], circles=var, add=T,...)}
    if (sum(y2,na.rm=T) > 0) {graphics::symbols(ordiscores[,1], ordiscores[,2], squares=y2, add=T,...)}
}

