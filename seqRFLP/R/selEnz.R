selEnz <-
function(names, enzdata = enzdata) {
    enz = enzdata[!is.na(match(enzdata[,1], names)), ]
    if( nrow(enz) == 0)
    { stop(paste("No data found in the enzdata named ", names, ".")) }
    return(enz)
}

