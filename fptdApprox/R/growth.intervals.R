growth.intervals <-
function (x, y, zeroSlope = 0.01) 
{
    l <- paste("0", paste(as.integer(diff(y) > 0), collapse = ""), 
        "0", sep = "")
    index.growth.from <- gregexpr("01", l)[[1]]
    if (index.growth.from[1] == -1) 
        return()
    index.growth.to <- gregexpr("10", l)[[1]]
    slope <- (180/pi) * atan((y[index.growth.to] - y[index.growth.from])/(x[index.growth.to] - 
        x[index.growth.from]))
    l <- slope > zeroSlope
    if (any(l))	    
        return(cbind(index.growth.from[l], index.growth.to[l]))
    else return()
}
