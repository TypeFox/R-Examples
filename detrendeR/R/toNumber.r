toNumber = function (x) {
    {
as.fractions = function (x){
as.integer((strsplit(x, "/")[[1]][1]))/ as.integer((strsplit(x, "/")[[1]][2]))
}
        if (length(grep("/", x)==1)) {
            X <- as.fractions(x)
        }
        else {
            X <- as.numeric(x)
        }
    }
    return(X)
}


#toNumber("2/3") 
#toNumber(0.44)
