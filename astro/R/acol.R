acol = function(col, alpha = 0){
    cols = col2rgb(col)/255
    return(rgb(t(rbind(col2rgb(col)/255)),alpha=alpha))
}

