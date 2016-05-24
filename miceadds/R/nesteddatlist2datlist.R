

nesteddatlist2datlist <- function(datlist){
    M <- length(datlist)
    B <- length(datlist[[1]])
    datlist0 <- as.list(1:(B*M))
    vv <- 1
    for (mm in 1:M){
        for (bb in 1:B){
            datlist0[[vv]] <- datlist[[mm]][[bb]]
            vv <- vv + 1
                        }
                    }
	class(datlist0) <- "datlist"
    return(datlist0)
            }