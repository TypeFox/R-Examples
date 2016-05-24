`count.diff.direction.error` <-
function(x,y){
    return(sum((diff(x) / diff(y))<0, na.rm=TRUE))
}

