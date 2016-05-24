objAndNames <- function(object, preferred, default)
{
  dimo <- dim(object)
  ndo <- length(dimo)
##
## 1.  NO 'dim' attribute
##
  if(ndo<1){
    n.o <- length(object)
    {
      if(is.list(preferred)){
        if((length(preferred)>0) && 
           (length(preferred[[1]])==n.o)){
          names(object) <- preferred[[1]]
          return(object)
        }
      }
      else if(length(preferred)==n.o){
        names(object) <- preferred
        return(object)
      }
    }
    {
      if(is.list(default)){
        if((length(default)>0) && 
           (length(default[[1]])==n.o) ){
          names(object) <- default[[1]]
          return(object)
        }
      }
      else if(length(default)==n.o){
        names(default) <- default
        return(object)
      }
    }
    return(object)
  }
##
## 2.  'dim' atribute
##
  dimn.o <- vector("list", ndo)
  for(id in 1:ndo){
    {
      if(is.list(preferred)){
        if((length(preferred)>=id) && 
           (length(preferred[[id]])==dimo[id])) {
          dimn.o[[id]] <- preferred[[id]]
          next
        }
      }
      else if((id==1) && (length(preferred)==dimo[id])){
        dimn.o[[id]] <- preferred
        next 
      }
    }
    {
      if(is.list(default)){
        if((length(default)>=id) && 
           (length(default[[id]])==dimo[id]) ){
          dimn.o[[id]] <- default[[id]]
          next 
        }
      }
      else if((id==1) && (length(default)==dimo[id]) ){
        dimn.o[[id]] <- default
        next 
      }
    }
  }
  dimnames(object) <- dimn.o
  object 
}


