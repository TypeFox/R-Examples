.removeZeroRows <-function(x){
    tmp <- rowSums(x)
    y = x[tmp > 0,, drop=FALSE ]
    return(y)
}
#' given matrix (columns protein rows peptides), compute minimal protein set using greedy algorithm
#' @param pepprot matrix as returned by prepareMatrix
#' @return list of protein peptide assignemnts (where winner takes all)
#' @export
#' @examples
#' library(prozor)
#'
#' data(protpepmetashort)
#' xx = prepareMatrix(protpepmetashort, weight= "count")
#' dim(xx)
#' es = greedy(as.matrix(xx))
#'
greedy <- function( pepprot ){
    ncolX = ncol(pepprot)
    res<-vector(ncolX , mode="list")
    idxx <-NULL
    for(i in 1:ncolX)
    {
        if(i %% 10 == 0){
            pepprot <- .removeZeroRows(pepprot)
            pepprot <- pepprot[,-idxx,drop=FALSE]
            idxx <-NULL
            message(paste("new dim" , dim(pepprot)))
        }
        if(nrow(pepprot) == 0){
            return(res[1:i])
        }
        oldtime <- Sys.time()
        pepsPerProt <- colSums(pepprot)
        if(max(pepsPerProt) == 0){
            return(res[1:i])
        }
        idx <- which.max(pepsPerProt)
        if(length(idx) > 1){
            idx<-idx[1]
        }
        idxx <- c(idxx, idx )
        dele <- pepprot[,idx]
        tmpRes = list(prot = colnames(pepprot)[idx], peps = rownames(pepprot)[dele > 0])
        res[[i]] <- tmpRes
        message(paste(i, " ", idx, " ", sum(dele)))
        if(sum(dele) > 0){
            set = cbind(rep(dele > 0, ncol(pepprot)))
            pepprot[set] <- 0
        }
        newtime <- Sys.time()
        message(paste("time ",newtime - oldtime))
    }
    return(res)
}
