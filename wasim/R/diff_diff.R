diff_diff<- function(x,y){
    return(sum((diff(x) / diff(y))<0, na.rm=TRUE))
}

