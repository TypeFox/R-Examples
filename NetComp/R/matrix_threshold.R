matrix_threshold<-function(matrix, threshold=NULL, minval=0, maxval=NULL, abs=TRUE, rmna=FALSE, ...){
    #x<-.Call("matrix_threshold",matrix, threshold, minval, maxval, abs, rmna, pkg="NetComp" )
    x<-.Call("matrix_threshold",matrix, threshold, minval, maxval, abs, rmna )
    if(is.null(row.names(matrix))==FALSE){
        row.names(x)<-row.names(matrix)
    }
     if(is.null(colnames(matrix))==FALSE){
        colnames(x)<-colnames(matrix)
    }
    x
}