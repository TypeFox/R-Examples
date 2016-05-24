mahaldis <- function(x) {

# check if numeric matrix and if no NA
if (!is.matrix(x) | any(is.na(x)) ) stop("x should be a matrix with no missing values")
x.labs <- rownames(x)

# pca function
pca <- function(Y, epsilon) {
   n <- nrow(Y)
   Y.svd = svd(Y)
   values = (1/(n-1))*Y.svd$d^2
   k <- sum(values > epsilon)
   values <- values[1:k]       
   U <- as.matrix(Y.svd$v[,1:k])
   F <- Y %*% U
   G <- sqrt(n-1)*as.matrix(Y.svd$u[,1:k])
   out <- list(values=values, U=U, F=F, G=G, k=k)
}

# center each column of x
x <- scale(x, center = TRUE, scale = FALSE)
n <- nrow(x)
epsilon <- sqrt(.Machine$double.eps)

# Mahalanobis distances among points are equal to Euclidean distances
# computed from the PCA matrix G, which positions the points in
# a PCA ordination space whose axes are stretched to account
# for the correlations among variables
G <- pca(x, epsilon)$G
D <- dist(G)
attr(D, "Labels") <- x.labs
attr(D, "Size") <- n
attr(D, "method") <- "Mahalanobis"
return(D)
}
