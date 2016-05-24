rjd <-
function(X, eps = 1e-06, maxiter = 100, na.action = na.fail)
    {
    X <- na.action(X)
    dim.X <- dim(X)
    
    if (length(dim.X)==2) type <- "Matrix"
    if (length(dim.X)==3) type <- "Array"
    if ((length(dim.X) %in% c(2,3))==FALSE) stop("'X' must have two or three dimensions")
    
    if (type == "Array")
        {
        if (dim.X[1] != dim.X[2]) stop("'X' must be an array with dim of the form c(p,p,k)")
        p <- dim.X[1]
        Xt <- aperm(X, c(1,3,2))
        X <- matrix(Xt, ncol=p)
        }
    
    if (!all(sapply(X, is.numeric)))  stop("'X' must be numeric")
    X <- as.matrix(X)
    X.data <- X
    
    p <- dim(X)[2]
    if (p==1) stop("'X' must be at least bivariate")
    kp <- dim(X)[1]
    
    
    k <- kp/p
    if (floor(k) != ceiling(k)) stop("'X' must me a matrix of k stacked pxp matrices")
    
    V <- diag(p)
    encore <- TRUE
    iter <- 0
    while (encore==TRUE)
        {
        iter <- iter +1
        encore <- FALSE
        for (i in 1:(p-1))
            {
            for (j in (i+1):(p))
                {
                Ii <- seq(i,kp,p)
                Ij <- seq(j,kp,p)
                
                g1 <- X[Ii,i]-X[Ij,j]
                g2 <- X[Ij,i]+X[Ii,j]
                
                g <- cbind(g1,g2)
                
                gg <- crossprod(g)
                ton <- gg[1,1]-gg[2,2]
                toff <- gg[1,2]+gg[2,1]
                theta <- 0.5*atan2(toff, ton+sqrt(ton*ton+toff*toff))
                
                cos.theta <- cos(theta)
                sin.theta <- sin(theta)
                
                if (abs(sin.theta)>eps)
                    {
                    encore <- TRUE
                    
                    Mi <- X[Ii,]
                    Mj <- X[Ij,]
                    
                    X[Ii,] <- cos.theta * Mi + sin.theta * Mj
                    X[Ij,] <- cos.theta * Mj - sin.theta * Mi
                    
                    col.i <- X[,i]
                    col.j <- X[,j]
                    
                    X[,i] <- cos.theta * col.i + sin.theta * col.j
                    X[,j] <- cos.theta * col.j - sin.theta * col.i
                    
                    temp <- V[i,]
                    V[i,] <- cos.theta * V[i,] + sin.theta * V[j,]
                    V[j,] <- cos.theta * V[j,] - sin.theta * temp
                    }
                }
            
            }
        if (iter >= maxiter) stop("maxiter reached without convergence")    
        }
    
        recomp <- function(X,V)
        {as.data.frame(V %*% tcrossprod(as.matrix(X), V))}
        Z <- split(as.data.frame(X.data), rep(1:k, each=p))
        Z2 <- sapply(Z, recomp, V=as.matrix(V), simplify=FALSE)
        Z3 <- matrix(unlist(lapply(Z2, t)), ncol=p, byrow=TRUE) 
        if (type == "Array")
            {
            D <- aperm(array(t(Z3), dim = c(p,p,k)), c(2,1,3), resize = FALSE)
            }
        else
            {
            D <- Z3
            }
        return(list(V=t(V),D=D))
    }
