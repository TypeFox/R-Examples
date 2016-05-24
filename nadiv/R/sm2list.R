#####################################
#adapted from code written by
#Jarrod Hadfield in the 
#MCMCglmm package
######################################
sm2list<-function(A, rownames = NULL, colnames=c("row", "column", "A"))
{
    ginv <- data.frame(Row = rep(1:length(A@p[-1]),
        diff(A@p)), Column = A@i + 1, Ainverse = A@x)
    ginv <- ginv[order(ginv$Row), ]
    ginv <- ginv[which(ginv$Row >= ginv$Column), ]
    attr(ginv, "rowNames") <- rownames
    names(ginv)<-colnames
    return(ginv)
}

