mpp <- function(data, gif, marks, params, gmap, mmap, TT){
    x <- list(data=data, gif=gif, marks=marks, params=params,
              gmap=gmap, mmap=mmap, TT=TT)
    class(x) <- "mpp"
    return(x)
}

