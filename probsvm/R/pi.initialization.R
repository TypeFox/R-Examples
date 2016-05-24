"pi.initialization" <- function(x, y, K, kernel, kparam, lambda, eps = 1.0e-10) {

    if (kernel == "radial") {
       kernel = "rbfdot"
       kpar = list("sigma" = kparam)
    } else if (kernel == "linear") {
       kernel = "vanilladot"
       kpar = list()
    } else if (kernel == "polynomial") {
       kernel = "polydot"
       kpar = list("degree" = kparam, "scale" = 1, "offset" = 1)
    } else stop("kernel function is not properly defined")

    cost <- 1/(2*lambda)
    temp <- ksvm(x = x, y = y, type = "C-svc", kernel = kernel, kpar = kpar, C = cost, tol = 1.0e-5, scaled = F)
    
    notRight <- temp@SVindex
    n<-length(y)
    alpha <- rep(0, n)
    alpha[notRight] <- unlist(temp@alpha) #* lambda
    
    Right <- unlist(setdiff(1:n, notRight))
    Left <- unlist(setdiff(notRight, which(alpha < (cost - 1.0e-10))))
    Elbow <- unlist(setdiff(notRight, Left))
    
    alpha <- alpha * lambda
    alpha0 <- temp@b * lambda
    
    gx <- t(alpha * y) %*% K[,Elbow]
    alpha0.elbow <- lambda * y[Elbow] - gx
    alpha0 <- mean(alpha0.elbow)

    obj <- list(alpha = alpha, alpha0 = alpha0, Elbow = Elbow, Left = Left, Right = Right)
   
obj
}
