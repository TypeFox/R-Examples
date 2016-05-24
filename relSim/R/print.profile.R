print.profile = function(x, horizontal = FALSE, ...){
    nLoci = length(x)/2
    i1 = 2*(1:nLoci) - 1
    i2 = i1 + 1
    strProf = sprintf("%2d/%-2d", x[i1], x[i2])

    if(horizontal){
        cat(paste(paste(strProf, collapse = ""), "\n"))
    }else{
        cat(paste(paste(strProf, collapse = "\n"), "\n"))
    }
}
