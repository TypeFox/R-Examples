#' Helper function to obtain a JSON-friendly format of a matrix formatted object
#' 
#' @param x A matrix.
#' @return A list for generating JSON format output.
#' @export

makeList<-function(x){
  if(ncol(x)>2){
    listSplit<-split(x[-1],x[1],drop=T)
    lapply(names(listSplit),function(y){list(name=y,children=makeList(listSplit[[y]]))})
  }else{
    lapply(seq(nrow(x[1])),function(y){list(name=x[,1][y],size=x[,2][y])})
  }
}