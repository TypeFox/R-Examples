strip = function(x, strip=" "){
    x = strsplit(x,"")[[1]]
    for(i in 1:length(strip)){
        s = strip[i]
        good = which(x!=s)
        if(length(good)>0){
            x = x[min(good):max(good)]
        }else{
            x = ""
        }
    }
    x = paste(x,collapse="")
    return(x)
}

