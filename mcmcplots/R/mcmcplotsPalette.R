mcmcplotsPalette <- function(n, type=c("rainbow", "sequential", "grayscale"), seq=NULL){
    type <- match.arg(type)
    if (type=="rainbow"){
        if(n==1)
            return(rainbow_hcl(1, start=240, l=50, c=100))
        return(rainbow_hcl(n, start=0, end=240, c=100))
    }
    if (type=="sequential"){
        return(sequential_hcl(n))
    }
    if (type=="grayscale"){
        return(gray((1:n/(n+1))))
    }
}
