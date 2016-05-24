`sample.median` <-
function(x){

    element1 <- c(NULL)
    element2 <- c(NULL)
    element3 <- c(NULL)


    replis <- unique(x[[4]][,"ID"])

    for (i in seq(along=replis)) {
        sample.lines <- which(x[[4]][,"ID"]==replis[i])
        sample.mtx <- x[[1]][sample.lines,]
        
    #### improvement 2010-10-13
            if (length(sample.lines) < 2) {

               print(paste("Warning: sample with identifier",replis[i],"has no replicate spots"))
               sample.mtx <- t(as.matrix(x[[1]][sample.lines,]))

             }
####

        medians <- apply(sample.mtx,2,median)
        element1 <- rbind(element1,medians)
        mads <- apply(sample.mtx,2,mad)
        element2 <- rbind(element2,mads)
        element3 <- rbind(element3,x[[4]][sample.lines[1],])
    }
    rownames (element1) <- c(1:nrow(element1))
    rownames (element2) <- c(1:nrow(element1))
    rownames (element3) <- c(1:nrow(element1))

    data.list <- list(expression=element1,
            error=element2,
            arraydescription=x[[3]],
            sampledescription=element3)


    return (data.list)
}

