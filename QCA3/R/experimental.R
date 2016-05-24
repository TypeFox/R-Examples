## BigCombnBySplit <- function(NR,k){
##     library(bigmemory)
##     if (NR < 2*k) stop("NR must be larger than 2*k.")
##     bigMatrix <- bigmemory:::filebacked.big.matrix(nrow=k,ncol=choose(NR,k),type="integer",backingfile="",descriptorfile="",init=0)
##     NR1 <-  seq(from=1,to=round(NR*0.5))
##     NR2 <-  seq(from=round(NR*0.5)+1,to=NR)

##     idx <- 1
##     for (i in 1:(k-1)){
##         combos1 <- combn(NR1, i)
##         combos2 <- combn(NR2, k-i)
##         for (j1 in 1:ncol(combos1)){
##             combos1n <- tryCatch(matrix(rep(combos1[,j1],each=ncol(combos2)),nrow=nrow(combos1),byrow=T),error=function(e)NULL)
##             if (!is.null(combos1n)){
##                 combos <- tryCatch(rbind(combos1n,combos2),error=function(e)NULL)
##                 if (!is.null(combos)){
##                     ## partially eliminate the for loop
##                     bigMatrix[,seq(from=idx,length=ncol(combos))] <- combos
##                     idx <- idx + ncol(combos)
##                 }} else {## partition combos2 and repeat the above will be better than loop ove all combos2
##                     ## or write to bigMatrix block by block directly
##                     for (j2 in 1:ncol(combos2)){
##                         bigMatrix[,idx] <- c(combos1[,j1],combos2[,j2])
##                         idx <- idx +1
##                     }
##                 }
##         }
##     }
##     if (length(NR1)>=k){
##         bigMatrix[,seq(from=idx,to=idx+(choose(length(NR1),k)-1))] <- combn(NR1, k)
##         idx <- idx+(choose(length(NR1),k))
##     }
##     if (length(NR2)>=k){
##         bigMatrix[,seq(from=idx,to=idx+(choose(length(NR2),k)-1))] <- combn(NR2, k)
##         ##idx <- idx+(choose(length(NR2),k))
##     }
##     bigMatrix
## }


## ## x=combnBySplit(100,5)
## ## combn(100,3) ~= combnBySplit(100,3)


## solveBigPIChart <- function (PIChart)
## {
##     require(biganalytics)
##     if (!is.logical(PIChart)) {
##         stop("Please use a logical matrix, such as an object returned by PIChart.\n")
##     }
##     if (all(dim(PIChart) > 1)) {
##         lpobj <- lpSolve:::lp(direction = "min", objective.in = rep(1,
##                                                  nrow(PIChart)), const.mat = t(PIChart), const.dir = ">=",
##                               1, all.bin = TRUE)
##         if (lpobj$status != 0)
##             stop("Can not solve this PMChart.")
##         k <- sum(lpobj$solution)
##         combos <- BigCombnBySplit(nrow(PIChart), k)
##         idx <- biganalytics:::apply(combos, 2, function(idx) all(colSums(PIChart[idx,,drop = FALSE]) > 0))
##         sol.matrix <- combos[,idx , drop = FALSE]
##     }
##     else {
##         sol.matrix <- matrix(seq_len(nrow(PIChart)), ncol = 1)
##     }
##     sol.matrix
## }
