acsh <- function (dotdata) {
    if (is.null(class(dotdata)) | class(dotdata) != "dotdata") {
        cat("Argument is not of class 'dotdata' \n")
        return(invisible())
    }
    nsp <- length(dotdata$Label)
    DistHomop <- matrix(0, nrow = nsp, ncol = nsp)
    if(nsp == 1) return(DistHomop) # Very trivial case  
    for (i in 1:(nsp-1)) {
      pt1 <- dotdata$occupancy[[i]]
      w1 <- dotdata$MSTsp[[i]]$wght
      for (j in (i+1):nsp) {
      		  pt2 <- dotdata$occupancy[[j]]
            w2 <- dotdata$MSTsp[[j]]$wght
            pts <- c(pt1, pt2[!(pt2 %in% pt1)])
            w <- w1
            w[pt1%in%pt2] <- w[pt1%in%pt2] + w2[pt2%in%pt1] 
            #This is valid because points are increasingly ordered by index number
            w <- c(w, w2[!(pt2 %in% pt1)])
            #Since weights by species points are normalized, the sum of vector w is equal 2. 
            DistHomop[i,j] <- DistHomop[j,i] <- 0.5*sum(w*apply(dotdata$PtSpDist[pts, c(i, j), drop = FALSE], 1, max))
          }
       }
    return(DistHomop)
}

