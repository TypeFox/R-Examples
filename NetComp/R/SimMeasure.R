SimMeasure<-function(data, threshold=NULL, ...){
    x<-.Call("SimMeasure",data, threshold, pkg="NetComp" )
    if(!is.null(row.names(data))){
        row.names(x)<-colnames(data); colnames(x)<-colnames(data)
    }
    x
}