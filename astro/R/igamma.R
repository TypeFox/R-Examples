igamma = function(x, s){
    # TRUE incomplete gamma function
    i = pgamma(x,s)*gamma(s)
    return(i)
}

