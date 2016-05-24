##------------------------
## Calculate h for multivariate method
##------------------------
checkh <- function(dat, h, fdelta, htype, centroids){
    n <- nrow(dat)
    p <- ncol(dat)

    if(centroids == "user"){
###################
## user
###################
        if(is.null(h)){
            if(fdelta == "mnorm"){
                if(htype == "AMISE"){
                    h <- AMISE(p, n)
                    cat("h missing. Using AMISE bandwidth.\n")
                }else if(htype == "ROT"){
                    h <- ROT(p, n)
                    cat("h missing. Using ROT bandwidth.\n")
                }else
                    stop("h is not given. Please provide a nonnegative value to h; OR set htype ='AMISE' or htype = 'ROT'.")
            }else{
                stop("h is not given. Please provide a nonnegative value to h")
            }
        }else{
            if(is.numeric(h) && h > 0)
                 htype <- "User specified"
            else{
                if(fdelta == "mnorm")
                    stop("h must be a nonnegative number; OR leave h = NULL, and set htype ='AMISE' or htype = 'ROT'.")
                else
                    stop("h must be a nonnegative number")
            }
        }
    }else if(centroids == "auto"){
###################
## auto
###################
        if(is.null(h)){
            if(fdelta == "mnorm"){
                if(htype == "AMISE"){
                    hstar <- AMISE(p, n)
                    h <- seq(hstar / 3, hstar * 3, length.out = 10)
                }
                else if(htype == "ROT"){
                    hstar <- ROT(p, n)
                    h <- seq(hstar / 3, hstar * 3, length.out = 10)
                    ## Fine
                }else{
                    stop("h is missing. htype is not 'AMISE' or 'ROT'.")
                }
            }else{
                stop("h is missing.")
            }
        }else{
            if(is.numeric(h) && all(h > 0))
                 htype <- "User specified"
            else{
                if(fdelta == "mnorm")
                    stop("h must be numberic or NULL. If h = NULL, htype must be one of 'AMISE' or 'ROT'. You can also set h = AMISE(p, n) or h = ROT(p, n), where p and n are dimension and number of observations.")
                else
                    stop("h must be a nonnegative number")
            }
        }
    }else
        stop("centroids only takes two options: 'user' or 'auto'.")
    
    ## if(!is.null(h)){
    ##     ##--------
    ##     ## h != NULL
    ##     ##--------    
    ##     if(h == "AMISE"){
    ##         if(fdelta == "mnorm"){
    ##             htype <- "AMISE"
    ##             h <- (4 / (p + 2)) ^ (1 / (p + 4)) * n ^ (-1 / (p + 4)) ## wand's AMISE
    ##         }else{
    ##             stop("h = 'AMISE' does not make sense if fdelta != 'mnorm'.")
    ##         }
    ##     }else if(h == "ROT"){
    ##         if(fdelta == "mnorm"){
    ##             htype <- "ROT"
    ##             h <- n ^ (-1 / (p + 4)) ## Scott's ROT
    ##         }else{
    ##             stop("h = 'ROT' does not make sense if fdelta != 'mnorm'.")
    ##         }
    ##     }else if(is.numeric(h) && h > 0){
    ##         ## good
    ##         htype <- "Specified"
    ##     }else{
    ##         stop("Please provide a nonnegative h value, or specify h = 'ROT' or h = 'AMISE'. If h is NULL, adpclust will attempt to find an h.")
    ##     }
    ## }else{
    ##     ##--------
    ##     ## h == NULL
    ##     ##--------
    ##     if(centroids == "user"){
    ##         if(fdelta == "mnorm"){
    ##             htype <- "AMISE"
    ##             h <- (4 / (p + 2)) ^ (1 / (p + 4)) * n ^ (-1 / (p + 4)) ## wand's AMISE
    ##             cat("\nBandwidth h not provided; Using AMISE value: h =", h, ".\n")
    ##         }else{
    ##             stop("For centroids == 'user', please provide a nonnegative h value. OR: if fdelta = 'mnorm', you can specify h = 'ROT' or h = 'AMISE'.")
    ##         }
    ##     }else{
    ##         ## fine
    ##     }
    ## }
    return(list(h = h, htype = htype))
}

##     ##------------------------
##     ## Define nclust in automatic variation if not given
##     ##------------------------
##     if(centroids == "auto" & is.null(nclust)) nclust = 2:10
    
##     ## This following function finds f and delta, given h
##     ## If h is missing then it is automatically found by highest silhouette.
##     if(centroids == "user" && fdelta == "mnorm"){
##         if(is.numeric(htype) && htype > 0)
##             temp.type <- "Given"
##         else if(h == "AMISE")
##             temp.type <- "AMISE"
##         else if(h == "ROT")
##             temp.type <- "ROT"
##         else
##             stop("'centroids is set to 'user'; Please provide a nonnegative h value, or specify h = 'ROT' or h = 'AMISE'.")
        
##         if(temp.type == "AMISE"){
##             h <- (4 / (p + 2)) ^ (1 / (p + 4)) * n ^ (-1 / (p + 4)) ## wand's AMISE
##         }else if(temp.type == "ROT"){
##             h <- n ^ (-1 / (p + 4)) ## Scott's ROT
##         }else{}
        
##         htype <- temp.type
##         fd.res <- findRd2(dat, h = h, plot = FALSE, dmethod = dmethod, htype = htype, 
##                           fdelta = fdelta, nclust = nclust)
##     }else if(centroids == "auto"){
        
##     }
## }
