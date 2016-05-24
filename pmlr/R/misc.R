tln<-function(x) log(ifelse(x<4.2e-307,4.2e-307,x))

texp<-function(x) exp(ifelse(x>709.7,709.7,x))

indicator.matrix <- function(x) {   
   n <- length(x)
   X <- matrix(data = 0, nrow = n, ncol = length(levels(x)) - 1, dimnames = list(NULL, levels(x)[-1])) 
   for (i in 2:length(levels(x))) X[,i - 1] <- x == levels(x)[i]
   X
}                        

vec <- function(x) cbind(as.vector(x))

# Horizontal direct product (i.e., *~ in GAUSS)

"%*~%" <- function (x, y) {
   m <- ncol(x)
   n <- ncol(y)
   x[, rep(1:m, each = n)] * x[, rep(1:n, m)]
} 
 