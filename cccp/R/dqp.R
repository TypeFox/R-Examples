##
## Function for creating an object of reference-class 'DQP'
dqp <- function(P, q, A = NULL, b = NULL, cList = list()){
   n <- ncol(P)
   if(is.null(A)){
       A <- matrix(0, nrow = 0, ncol = n)
    } 
    if(is.null(dim(A))){
        A <- matrix(A, nrow = 1)
    }
    if(is.null(b)){
        b <- numeric()
    } 
    if(length(cList) > 0){
        cone <- unlist(lapply(cList, function(x) x[["conType"]]))
        if(!all(cone %in% c("NNOC", "SOCC", "PSDC"))){
            stop("List elements of cone constraints must be either created by calls to:\n'nnoc()', or 'socc()', or 'psdc()'.\n")
        }
        K <- length(cList)
        GList <- lapply(cList, function(x) x[["G"]])
        G <- do.call("rbind", GList)
        h <- do.call("rbind", lapply(cList, function(x) x[["h"]]))
        ridx <- cumsum(unlist(lapply(GList, nrow)))
        sidx <- cbind(c(0, ridx[-length(ridx)]), ridx - 1)
        dims <- as.integer(unlist(lapply(cList, function(x) x[["dims"]])))
        cList <- new(CONEC, cone, G, h, sidx, dims, K, n)
    } else {
        cList <- new(CONEC, as.integer(n))
    }
   ans <- new(DQP,
              P = P,
              q = q,
              A = A,
              b = b,
              cList = cList)
   return(ans)
}
