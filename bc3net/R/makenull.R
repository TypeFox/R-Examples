

makenull <- function(dataset,nullit=NA, estimator="pearson", disc="equalwidth"){
     
     if(!is.matrix(dataset)){
        cat("dataset not a matrix\n")
        dataset <- as.matrix(dataset)
     }
     
     rown <- nrow(dataset)
     coln <- ncol(dataset)
     null <- c()

     if(is.na(nullit)){
       nullit=ceiling(10^5/(((rown*rown)/2)-rown))
     }     

     for (i in 1:nullit) {
         expdata <- matrix(sample(dataset), rown, coln)
         mim <- mimwrap(expdata, estimator=estimator, disc=disc)
         miv <- mim[upper.tri(mim)]
         miv <- miv[miv != 0]
         miv <- as.vector(miv)
         null <- c(null, miv)
     }

     null <- sort(null)
     return(null)
}
