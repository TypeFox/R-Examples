myFrac <-
function(num, denom = NULL){
  if(is.null(denom)){
    temp <- strsplit(as.character(num),"/")
    num <- lapply(temp,function(x){x[1]})
    denom <- lapply(temp,function(x){x[2]})
  }
  paste0("^",num,"^/~",denom,"~")
}
