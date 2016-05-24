smoother <- function(notes, method ="median", order=4, times=2){
    require("pastecs")
    notes[is.na(notes)] <- 999
    dmnotes <- unclass(pastecs::decmedian(notes, order = order, times = times)$series[,1])
    is.na(dmnotes[dmnotes==999]) <- TRUE
    return(as.numeric(dmnotes))
}
