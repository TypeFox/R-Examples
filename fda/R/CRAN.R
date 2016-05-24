CRAN <- function(CRAN_pattern, n_R_CHECK4CRAN){
##
## 1.  get environment variables
##
    gete <- Sys.getenv()
    ngete <- names(gete)
    i <- seq(along=gete)
##
## 2.  check CRAN_pattern
##
    if(missing(CRAN_pattern)){
        if('_CRAN_pattern_' %in% ngete){
            CRAN_pattern <- gete['_CRAN_pattern_']
        } else CRAN_pattern <- '^_R_'
    }
##
## 3.  check n_R_CHECK4CRAN
##
    if(missing(n_R_CHECK4CRAN)){
        if('_n_R_CHECK4CRAN_' %in% ngete){
            n_R_CHECK4CRAN <- as.numeric(gete['_n_R_CHECK4CRAN_'])
        } else n_R_CHECK4CRAN <- 5
    }
##
## 4.  Check
##
    for(pati in CRAN_pattern)
        i <- i[grep(pati, ngete[i])]
##
## 5.  Done
##
    cran. <- (length(i) >= n_R_CHECK4CRAN)
    attr(cran., 'Sys.getenv') <- gete
    attr(cran., 'matches') <- i
    cran.
}


