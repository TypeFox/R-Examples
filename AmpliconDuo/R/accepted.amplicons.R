accepted.amplicons <-
function(x){
  
  #### internal function #######
  get.accepted.amplicons <-
    function(x){
      accep <- which(x$rejected == FALSE)
      return(as.integer(row.names(x[accep,])))
    }
  ######################
  if(class(x) == "list"){
    acamp <- lapply(X = x, FUN = get.accepted.amplicons)
    names(acamp) <- names(x)
  }else{
    acamp <- get.accepted.amplicons(x)
  }
  return(acamp)
}
