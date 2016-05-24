BASIX.combnapply <- function(vec,mode="*"){

ret <- .Call("combnapply_C",vec,mode)

if(mode=="=="){
ret <- as.logical(ret)
}
return(ret)

}
