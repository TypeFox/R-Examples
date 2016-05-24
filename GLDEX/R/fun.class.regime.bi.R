"fun.class.regime.bi"<-
function (data, perc.cross, fun.cross) 
{
 
    if (is.function(fun.cross)) {
        data <- sort(data)
        index <- fun.cross(data, 2)$clustering
        data.a <- sort(data[index == 1])
        data.b <- sort(data[index == 2])
    }
    if (is.logical(fun.cross)) {
        data.a <- sort(data[fun.cross == TRUE])
        data.b <- sort(data[fun.cross == FALSE])
    }
    if (perc.cross != 0) {
        if (ceiling(perc.cross * length(data.a)) > length(data.b) || 
            ceiling(perc.cross * length(data.b)) > length(data.a)) {
            stop("Please decrease perc.cross level, current level is too high")
        }
        data.a <- c(data.a, data.b[sample(1:length(data.b), ceiling(perc.cross * 
            length(data.a)) - 1)])
        data.b <- c(data.b, data.a[sample(1:length(data.a), ceiling(perc.cross * 
            length(data.b)) - 1)])
        if (is.element(min(data), data.a) == FALSE) {
            data.a <- c(data.a, min(data))
        }
        if (is.element(max(data), data.a) == FALSE) {
            data.a <- c(data.a, max(data))
        }
        if (is.element(min(data), data.b) == FALSE) {
            data.b <- c(data.b, min(data))
        }
        if (is.element(max(data), data.b) == FALSE) {
            data.b <- c(data.b, max(data))
        }
    }
    return(list("data.a"=data.a, "data.b"=data.b))
}



