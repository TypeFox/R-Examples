ESdist <-
function(adate)

{
    # estimate Earth-Sun distance for date adate
    # in format "YYYY-MM-DD"
    # result is in AU

    edist <- julian(as.Date(adate), origin=as.Date(paste(substring(adate, 1, 4), "12", "31", sep="-")))[[1]]
    edist <- 1 - 0.016729 * cos((2*pi) * (0.9856 * (edist - 4)/360))
    
    edist
}

