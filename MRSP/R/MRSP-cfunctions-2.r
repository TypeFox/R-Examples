################################################################################
#####    C functions for MRSP.                                             #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 28.02.2014, 20:32                                  #####
################################################################################

##### a more efficient function for computing the rowwise maximum of a matrix
maxRow <- cmpfun(function(x){
 if (is.data.frame(x)) x <- as.matrix(x)
 if (!is.array(x) || length(dn <- dim(x)) < 2)
      stop("`x' must be an array of at least two dimensions")
 p <- ncol(x)
 dn <- nrow(x)
 out <- dn
 if (!is.double(x)) x <- as.double(x)
 .C("maxRow", x, as.integer(dn), as.integer(p), double(out),
     DUP = TRUE, PACKAGE = "MRSP")[[4]]
})


##### an efficient function for computing the cumsum of each row in a matrix
rowCumsum <- cmpfun(function(x){
 if (is.data.frame(x)) x <- as.matrix(x)
 if (!is.array(x) || length(dim(x)) < 2)
      stop("`x' must be an array of at least two dimensions")
 p <- ncol(x)
 n <- nrow(x)
 if (!is.double(x)) x <- as.double(x)
 out <- x
 matrix(.C("rowCumsum", x, as.integer(n), as.integer(p), out,
            DUP = TRUE, PACKAGE = "MRSP")[[4]], nrow = n)
})

##### an efficient function for computing the cumprod of each row in a matrix
rowCumprod <- cmpfun(function(x){
 if (is.data.frame(x)) x <- as.matrix(x)
 if (!is.array(x) || length(dim(x)) < 2)
      stop("`x' must be an array of at least two dimensions")
 p <- ncol(x)
 n <- nrow(x)
 if (!is.double(x)) x <- as.double(x)
 out <- x
 matrix(.C("rowCumprod", x, as.integer(n), as.integer(p), out,
            DUP = TRUE, PACKAGE = "MRSP")[[4]], nrow = n)
})

##### a function for constructing the modified response matrix in gradient
yseqlog.constructor <- cmpfun(function(xn, xd, m){
 if (is.data.frame(xn)) x <- as.matrix(xn)
 if (!is.array(xn) || length(dim(xn)) < 2)
      stop("`xn' must be an array of at least two dimensions")
 p <- ncol(xn)
 n <- nrow(xn)
 if(any(dim(xn) != dim(xd))) stop("xn and xd do not match in yseqlogconstructor")
 if (!is.double(xn)) x <- as.double(xn)
 if (!is.double(xd)) x <- as.double(xd)
 out <- xn
 matrix(.C("yseqlogconstructor", xn, xd, as.integer(n), as.integer(p), as.integer(m), out,
            DUP = TRUE, PACKAGE = "MRSP")[[6]], nrow = n)
})


##### a function for computing the lasso as used in MRSP, that is, a rowwise L2-norm shrinkage that
##### reduces to the ordinary lasso if predictors are not grouped by grpindex. 
lassoC <- cmpfun(function(u, lambda, w, dp, dn, norms, out){
.C("lassoC", u, lambda, w, as.integer(dp), as.integer(dn), norms, out, DUP = TRUE, PACKAGE = "MRSP")
})

