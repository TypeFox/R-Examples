gdist.total <- function (longlat, units = 'nm', segments = TRUE, digits = 2) 
{
    if (!(is.data.frame(longlat) | is.matrix(longlat)) | nrow(longlat) < 2) 
        return(0)

    SEGMENTS <- apply(cbind(longlat, rbind(longlat[2:nrow(longlat), ], c(NA, NA))), 
        1, function(x) {
               gdist(x[1], x[2], x[3], x[4], units = units)
        })
    N <- length(SEGMENTS)

    if(segments & N > 2) {
       
       print(data.frame(Segment = 1:(N-1), Distance = paste(round(SEGMENTS[1:(N-1)], digits), units)))
       flush.console()
    }

    sum(SEGMENTS, na.rm = TRUE)

}


