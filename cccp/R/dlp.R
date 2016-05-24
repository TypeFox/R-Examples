##
## Function for creating an object of reference-class 'DLP'
dlp <- function(q, A = NULL, b = NULL, cList = list()){
    if(is.matrix(q)){
        warning("Matrix provided for q, extracting first column for argument 'q'.\n")
        q <- q[, 1]
    }
    n <- length(q)
    if(is.null(A)){
        A <- matrix(0, nrow = 0, ncol = n)
    } 
    if(is.null(dim(A))){
        A <- matrix(A, nrow = 1)
    }
    if(is.null(b)){
        b <- numeric()
    }
    K <- length(cList)
    if(K < 1){
       warning("LP in standard form: Adding non-negativity constraint(s).\n")
       G <- -diag(n)
       h <- rep(0, n)
       cList <-  list(nnoc(G = G, h = h))
       K <- 1
   }
   if(K > 0){
       cone <- unlist(lapply(cList, function(x) x[["conType"]]))
       if(!all(cone %in% c("NNOC", "SOCC", "PSDC"))){
           stop("List elements of cone constraints must be either created by calls to:\n'nnoc()', or 'socc()', or 'psdc()'.\n")
       }
       GList <- lapply(cList, function(x) x[["G"]])
       G <- do.call("rbind", GList)
       h <- do.call("rbind", lapply(cList, function(x) x[["h"]]))
       ridx <- cumsum(unlist(lapply(GList, nrow)))
       sidx <- cbind(c(0, ridx[-length(ridx)]), ridx - 1)
       dims <- as.integer(unlist(lapply(cList, function(x) x[["dims"]])))
       cList <- new(CONEC, cone, G, h, sidx, dims, K, n)
    } else {
        stop("LP only with equality constraints; undefined or exact solution.\n")
    }
    ans <- new(DLP,
               q = q,
               A = A,
               b = b,
               cList = cList)
    return(ans)
}
