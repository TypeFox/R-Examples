#projection matrix
projMat <-
function(X) X %*% ginv(t(X) %*% X) %*% t(X)
