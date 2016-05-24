"process.args" <- function(...){
  a <- list(...)
  attr <- attributes(a[[which.max(unlist(lapply(a,length)))]])
  a <- lapply(a,as.vector)
  out <- do.call("rbind",a)
  
  out <- split(out,row(out))
  names(out) <- paste("arg",1:length(a),sep="")
  return(c(out,attr=list(attr)))
}

#"process.2.args" <- function(a1,a2){
#  lens <- c(length(a1),length(a2))
#  attributes.list <- list(attributes(a1),attributes(a2))
#  attributes.wanted <- attributes.list[[which.max(lens)]]
#
#  jj <- rbind(as.vector(a1),as.vector(a2))
#  return(list(arg1=as.vector(jj[1,]),
#              arg2=as.vector(jj[2,]),
#              attr=attributes.wanted
#              )
#         )
#}
#
#process.3.args <- function(a1,a2,a3){
#  lens <- c(length(a1),length(a2),length(a3))
#  attributes.list <- list(attributes(a1),attributes(a2),attributes(a3))
#  attributes.wanted <- attributes.list[[which.max(lens)]]
#
#  jj <- rbind(as.vector(a1),as.vector(a2),as.vector(a3))
#  return(list(arg1=as.vector(jj[1,]),
#              arg2=as.vector(jj[2,]),
#              arg3=as.vector(jj[3,]),
#              attr=attributes.wanted
#              )
#         )
#}
#process.4.args <- function(a1,a2,a3,a4){
#  lens <- c(length(a1),length(a2),length(a3),length(a4))
#  attributes.list <- list(attributes(a1),attributes(a2),attributes(a3),attributes(a4))
#  attributes.wanted <- attributes.list[[which.max(lens)]]
#
#  jj <- rbind(as.vector(a1),as.vector(a2),as.vector(a3),as.vector(a4))
#  return(list(arg1=as.vector(jj[1,]),
#              arg2=as.vector(jj[2,]),
#              arg3=as.vector(jj[3,]),
#              arg4=as.vector(jj[4,]),
#              attr=attributes.wanted
#              )
#         )
#}
#
#
strictify <- function(val,status)
  {
    val[status>0] <- NaN
    return(val)
  }
