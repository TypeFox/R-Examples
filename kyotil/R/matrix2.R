# matrix functions used to improve performance of krm.score.test
# originally from Krisztian Sebestyen ksebestyen@gmail.com

.as.double <- function(x, stripAttributes=FALSE){
    if(!stripAttributes){
        if(is.double(x)) return(x)   # no duplicate copy
        storage.mode(x) <- 'double'  # yes duplicate copy
        return(x)
    }
    return(as.double(x))             # duplicate copy if has attributes or not double 
}

# t(X) %*% diag(d) %*% X
# X is a nxp matrix, D is a diagonal matrix of n by n or an array of length n
tXDX <- function(X,D){ 
    
    # reduced D to an array
    if(NROW(D) == NCOL(D)) D <- diag(D)     
    n<-length(D)  
    stopifnot(NROW(X) == n)
    
    if(!is.matrix(X)) X <- as.matrix(X)
        
    D <- .as.double(D) # critical to use .as.double here, otherwise it will try to copy
    X <- .as.double(X) # critical to use .as.double here, otherwise it will try to copy
    .Call("xdx", X, D)
}


DXD <- function(d1,X,d2){ 
    
    if(NROW(d1) == NCOL(d1)) d1 <- diag(d1) # reduced to diagonal of the matrix
    if(NROW(d2) == NCOL(d2)) d2 <- diag(d2) # reduced to diagonal of the matrix
    
    n<-length(d1)
    
    if((NROW(X) != n) || (NCOL(X) != n)) 
    if(!is.matrix(X)) X <- as.matrix(X)
        
    # these are needed b/c dxd2 expects input type of double (by using REAL) and do not perform type coersion
    d1 <- .as.double(d1) # critical to use .as.double here, otherwise it will try to copy
    d2 <- .as.double(d2)
    X  <- .as.double(X)
    
#    .C("dxd",n,d1,X,d2,dxd = X*0.0,DUP = FALSE,NAOK = FALSE)$dxd
    
    .Call("Call_dxd", d1, X, d2)
}


symprod <- function(S,X){

    if(!is.matrix(S) || !is.matrix(X)) stop("Both 'S' and 'X' have to be matrices.")
    S <- .as.double(S)
    X <- .as.double(X)
    #Y <- X * 0.0
    if(
        (NCOL(S) != NROW(X)) || 
        (NCOL(S) != NROW(S))  
        #|| !all(dim(X) == dim(Y))
    ) stop("Dimension mismatch")
    
#    res <- .C( 
#        "symprod",  
#        M = NROW(Y),
#        N = NCOL(Y),
#        A = .as.double(S),
#        B = as.double(X),
#        C = Y
#    ) 
#
#    # res <- .Fortran( 
#        # "dsymm",  
#        # SIDE = 'L',
#        # UPLO = 'U',
#        # M = NROW(Y),
#        # N = NCOL(Y),
#        # alpha = as.double(1),
#        # A = .as.double(S),
#        # LDA = NROW(S),
#        # B = as.double(X),
#        # LDB = NROW(X),
#        # beta = as.double(0),
#        # C = Y,
#        # LDC = NROW(Y) 
#    # )     
#    res$C
    
    .Call("Call_symprod", S, X)

}


txSy <- function(x,S,y){
    n <- length(x)
    
    if( (n != NROW(S)) || (n != NCOL(S)) || (n != length(y)) )stop("Dimension mismatch")
    
    S <- .as.double(S)
    x <- .as.double(x)
    y <- .as.double(y)


    #.C('txSy',length(x),.as.double(S),.as.double(x),.as.double(y),double(n),out=double(1))$out
    .Call('Call_txSy', x, S, y)
}

# rep matrix along rows in C
# times is an integer >= 0
# each is >= 0, an integer or a vector of length nrow(x)
rrbind <- function(x,times = 1,each = 1){
    if(!is.matrix(x)) stop("'x' must be a matrix")
    vec.each <- NULL
    if(!length(each)) each <- 0
    each <- as.integer(each)
    if(length(each) > 1) vec.each <- rep(each,length.out = nrow(x)) 
    if(!length(times)) times <- 0   
    .Call('Call_rrbind',.as.double(x),as.integer(times),each,vec.each)  
}

# rep matrix along columns in C
# times is an integer >= 0
# each is >= 0, an integer or a vector of length ncol(x)
rcbind <- function(x,times = 1,each = 1){
    if(!is.matrix(x)) stop("'x' must be a matrix")
    vec.each <- NULL
    if(!length(each)) each <- 0
    each <- as.integer(each)
    if(length(each) > 1)vec.each <- rep(each,length.out = ncol(x))  
    if(!length(times)) times <- 0   
    .Call('Call_rcbind',.as.double(x),as.integer(times),each,vec.each)  
}

# x: {m x n} matrix
# k: integer in {-m , .. , 0 , .. , n} determining position of diagonal 
# k = 0   <-> x[1,1] = diag[1]
# k = 1   <-> x[1,2] = diag[1]
# k = -1  <-> x[2,1] = diag[1]
# diag: what the 'diagonal' of 'x' should be set to, 0 if NULL

lower.trap <- function(x,pos = 0,diag = NULL){
    if(!is.matrix(x)) x <- as.matrix(x)
    if(!is.null(diag))
        diag <- rep(.as.double(diag),length.out = min(dim(x)))
    .Call('Call_lower_trap',.as.double(x),diag,as.integer(pos)) 
}

upper.trap <- function(x,pos = 0,diag = NULL){
    if(!is.matrix(x)) x <- as.matrix(x)
    if(!is.null(diag))
        diag <- rep(.as.double(diag),length.out = min(dim(x)))
    .Call('Call_upper_trap',.as.double(x),diag,as.integer(pos)) 
}
