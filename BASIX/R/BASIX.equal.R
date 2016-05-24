BASIX.equal <- function(a,b){


if(is.character(a)){
eq <- .Call("Ccompare2",a,b,PACKAGE="BASIX")
}else{
eq <- .Call("Ccompare",a,b,PACKAGE="BASIX")
}

eq <- as.logical(eq)
return(eq)

}
