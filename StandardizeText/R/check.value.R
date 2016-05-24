check.value <-
function(value,type) {
  if(type!="boolean" && type!="string" && type!="integer" && type!="str_int" && type!="vector") {
    return()
  }
  if(is.vector(value)) {
    if(type=="vector" && !is.list(value)) {
      return(TRUE)
    } else if(length(value)!=1) {
      return(FALSE)
    }else if (type=="boolean" && is.logical(value)) {
      return(TRUE)
    }else if((type=="string" || type=="str_int") && is.character(value)) {
      return(TRUE)
    }else if((type=="integer" || type=="str_int") && is.numeric(value) && as.integer(value)==value) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else if(type=="vector" && is.factor(value)) {
    return(TRUE)
  }else if(is.null(value)) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}
