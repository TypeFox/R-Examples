# the function polyS() is like the poly() S-plus function
# it works with conjuction to the poly.matrix() function and
# allows to create orthigonal polynomials in more that one dimension

polyS<- function(x, ...) # like poly() in S
{
    vars <- list(...)   
    #call with single arguments, no degree i.e. polyS(x)
    if(length(vars) == 0) return(poly.matrix(x, degree = 1))    
    #call with single arguments and degree equal to value  i.e. polyS(x,2)
    if(length(vars) == 1 & length(vars[[1]]) == 1)
        return(poly.matrix(x, degree = vars[[1]]))  #  
    #call with 2 or more variables i.e polyS(x1,x2,x3,2)
    mat <- Recall(...)
    ad <- attr(mat, "degree")
    degree <- max(ad)
    mat0 <- poly.matrix(x, degree)
    ad0 <- attr(mat0, "degree")    
     ncol.mat <- if(!is.matrix(mat))  1 else dim(mat)[[2]]
    ncol.mat0 <- if(!is.matrix(mat0)) 1 else dim(mat0)[[2]]
    index0 <- rep(1:ncol.mat0, ncol.mat)
    index <- rep(1:ncol.mat, rep(ncol.mat0, ncol.mat))
    newad <- ad0[index0] + ad[index]
    retain <- newad <= degree
    if(any(retain))
        newmat <- mat0[, index0[retain], drop = FALSE] * mat[, index[retain], drop = FALSE]
    else newmat <- NULL
    d0 <- dimnames(mat0)
    d <- dimnames(mat)
    dn <- d0
    dn[[2]] <- paste(d0[[2]][index0[retain]], d[[2]][index[retain]], sep = ".")
    dimnames(mat0) <- d0
    dimnames(mat) <- d
    if(length(newmat))
        dimnames(newmat) <- dn
    newmat <- cbind(mat0, mat, newmat)
    newad <- c(ad0, ad, newad[retain])
    attr(newmat, "degree") <- newad
    newmat
}




poly.matrix <- function(m, degree = 1) # like poly.matrix in S
{
    dd <- dim(m)
    if(is.null(dd))
         poly(m, degree) #test
    else {
        class(m) <- NULL
         if(is.list(m)) {
            m <- as.vector(m)
            names(m) <- NULL
            do.call("polyS", c(m, list(degree)))
           }
            else {
            if((dd2 <- dd[2]) == 1)
                poly(m, degree)
            else {
                listm <- vector("list", dd2)
                for(i in 1:dd2)
                  listm[[i]] <- m[, i]
                do.call("polyS", c(listm, degree))
                 }
           }
        }
}
