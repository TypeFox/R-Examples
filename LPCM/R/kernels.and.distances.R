

# computes all distances between a set of vectors and another vector.
# (up to the constant sqrt(d)) 
 distancevector <- function(X, y, d = "euclid", na.rm = TRUE) {

      X <- as.matrix(X)
          if (dim(X)[2]==1 && dim(X)[1]==length(y)){X<-t(X)}
       y <- as.numeric(y)
   
        vdisseuclid <- function(X, y, na.rm = TRUE) { # from R package hopach.
           if (!is.matrix(X)) 
              stop("First arg to vdisseuclid() must be a matrix")
          if (!is.vector(y)) 
             stop("Second arg to vdisseuclid() must be a vector")
          dX <- dim(X)
          p <- dX[1]
          n <- dX[2]
          if (length(y) != n) 
              stop("Matrix and vector dimensions do not agree in vdisseuclid()")
          if (na.rm) {
             N1 <- rowSums(!is.na(X))
             N2 <- sum(!is.na(y))
             N3 <- (!is.na(X)) %*% (!is.na(y))
             X[is.na(X)] <- 0
             y[is.na(y)] <- 0
             suppressWarnings(out <- sqrt(as.vector(rowSums(X^2)/N1 + 
                 sum(y^2)/N2 - 2 * X %*% y/N3)))
           }
           else suppressWarnings(out <- sqrt(as.vector(rowMeans(X^2) + 
               mean(y^2) - 2 * X %*% y/n)))
           out[out == "NaN"] <- 0
           return(out)
         }
        
        if (d == "euclid") {
           return(vdisseuclid(X, y, na.rm))
         }
        stop("Distance metric ", d, " not available")
    }


# computes row by row the pairwise distances between two sets of vectors.
vecdist <- function(X,Y){
      if (is.data.frame(X)){
           X<- as.matrix(X)
      } else if (!is.matrix(X)){
           X<-matrix(X, nrow=1)
      }
      if (is.data.frame(Y)){
           Y<- as.matrix(Y)
      } else if (!is.matrix(Y)){
        Y<-matrix(Y, nrow=1)
      }
      sqrt( apply(  (X-Y)^2 , 1,sum))     
 }


# computes the minimal distance between a vector and a set of vectors 
 mindist <- function(X, y){
        d <- length(y)
        s <- sqrt(d) * distancevector(X, y)
        return(list(mindist = min(s), closest.item = order(s)[1]))
    }


# Univariate kernel
kern <- function(y, x = 0, h = 1){
    1/h * dnorm((x - y)/h)
}

# Multivariate kernel 
kernd <- function(X,x,h){
   if (!is.matrix(X) && !is.data.frame(X)){ X<- matrix(X, nrow=1)}
   x<-as.numeric(x)
   d<-length(x)
   k<-1
   for (j in 1:d){k<- k* kern(X[,j],x[j],h[j])}
   k
   }

# Euclidian norm
enorm <- function (x){
         sum(x^2)
     }

# Multivariate kernel density estimator

kdex <-  function(X, x, h){
          n<-dim(X)[1]
          1/n*sum(kernd(X,x,h))
        }
