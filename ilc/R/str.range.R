str.range <-
function(x){
    if (is.numeric(x)) x <- sort(x)
    paste(x[1], x[length(x)], sep=' - ')
}
