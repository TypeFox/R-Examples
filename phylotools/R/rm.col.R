#### Function rm.col as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


rm.col <-
function(mat, to.rm = NULL){
    if(is.null(to.rm)){
        return(mat)
    } 
    else{
        if(!is.character(to.rm)){
            stop("The column names must be in class character.")
        }
        cnames <- colnames(mat)
        if(!all(to.rm%in%cnames)){
           notinclude <- to.rm[which(!to.rm%in%cnames)]
           cat("Column", notinclude, "\n", 
           "not found in column names of the mat, ignored\n")
        }
        ind <- -which(cnames%in%to.rm)
        res <- mat[,ind]
    }
    return(res)
}

