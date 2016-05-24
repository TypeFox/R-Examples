deg2dec <-
function(h,m,s){
    if(h < 0){
        m = - m
        s = -s
    }
    res = h + m/60 + s/3600
    return(res)
}

