LP.comoment <- function(x,y,zero.order = TRUE,m = 5){
    if(length(x) != length(y))
        stop("\n Different lengths of x and y \n")
    Tx <- LPTrans(x  ,m)
    Ty <- LPTrans(y  ,m)
    lx <- cor(x,Ty)
    ly <- cor(y,Tx)
    com.matrix <- cbind(c(1,lx), rbind(ly,cov(Tx, Ty)))
    rownames(com.matrix) <- paste("T.",0:m,"(X)", sep ="")
    colnames(com.matrix) <- paste("T.",0:m,"(Y)", sep = "")
    if(zero.order == FALSE)
        com.matrix <- com.matrix[,-1][-1,]
    return(com.matrix)
    }
