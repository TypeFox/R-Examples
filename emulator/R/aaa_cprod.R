ht <- function(x){  #Hermitian transpose
 if(is.complex(x)){
    return(t(Conj(x)))
 } else {
    return(t(x))
 }
}

"cprod" <- function(x,y=NULL){
 if(is.complex(x) | is.complex(y)){
    if(is.null(y)){
     return(crossprod(Conj(x),x))
    } else {
     return(crossprod(Conj(x),y))
    }
 } else {
    return(crossprod(x,y))
 }
}

"tcprod" <- function(x,y=NULL){
 if(is.complex(x) | is.complex(y)){
    if(is.null(y)){
     return(tcrossprod(x,Conj(x)))
    } else {
     return(tcrossprod(x,Conj(y)))
    }
 } else {
    return(tcrossprod(x,y))
 }
}

"quad.form" <-
function (M, x, chol = FALSE)
{
    if (chol == FALSE) {
       return(drop(cprod(cprod(M, x), x)))
    }
    else {
       jj <- cprod(M, x)
       return(drop(cprod(jj, jj)))
    }
}

"quad.form.inv" <-
function (M, x)
{
    drop(cprod(x, solve(M, x)))
}

"quad.3form" <-
function(M,left,right)
{
 cprod(cprod(M,left),right)
}  

"quad.3tform" <- function(M,left,right)
{
  tcprod(left,tcprod(right,M))
}

"quad.tform" <-
function(M,x)
{  
 tcprod(tcprod(x,M),x)
}

"quad.tform.inv" <-
function(M,x){
 drop(quad.form.inv(M,ht(x)))
}

"quad.diag" <- 
function(M,x){   # algebraically, diag(quad.form(M,x))
    colSums( cprod(M,x) * Conj(x))
}

"quad.tdiag" <-
function(M,x){   # algebraically, diag(quad.tform(M,x))
    rowSums( tcprod(x,M) * Conj(x))
}

#"cmahal" <- 
#    function (z, center, cov, inverted = FALSE, ...) 
#{
#    if(is.vector(z)){
#        z <- matrix(z, ncol = length(z))
#    } else {
#        z <- as.matrix(x)
#    }
#    
#    if (!inverted) { cov <- solve(cov, ...)}
#    quad.diag(cov,sweep(z, 2, center))
#}

