cjd <- function(X, eps=1e-06, maxiter = 100){

    dim.X<-dim(X)

    if (length(dim.X)==2) type<-"Matrix"
    if (length(dim.X)==3) type<-"Array"
    if ((length(dim.X) %in% c(2,3))==FALSE) stop("'X' must have two or three dimensions")

    if (type == "Matrix")
     {
      p <- dim.X[2]
      K <- dim.X[1]/p
      if (floor(K) != ceiling(K)) stop("'X' must be a matrix of k stacked pxp matrices")
      X <- array(t(X), c(p, p, K))
      dim.X <- dim(X)
     }


    if (dim.X[1] != dim.X[2]) stop("'X' must be an array with dim of the form c(p, p, K)")
    p <- dim(X)[1]
    k <- dim(X)[3]

    kp <- k*p
    Xt <- aperm(X, c(1, 3, 2))
    X <- matrix(Xt, nrow = p, byrow=TRUE)

    b <- complex(real=c(1, 0, 0, 0, 1, 1, 0, 0, 0), imaginary=c(0, 0, 0, 0, 0, 0, 0, -1, 1))
    B <- matrix(b, 3, 3, byrow=TRUE)
    Bt <- Conj(t(B))

    # initial value
    V <- diag(p)
    encore <- TRUE
    iter <- 0

    while (encore == TRUE) {
        iter <- iter + 1
        encore <- FALSE

        for (i in 1:(p - 1)) {
            for (j in (i + 1):(p)) {
                Ii <- seq(i, kp, p)
                Ij <- seq(j, kp, p)
                g <- rbind(X[i, Ii] - X[j, Ij], X[i, Ij], X[j, Ii])
                gg <- B %*% g %*% Conj(t(g)) %*% Bt
                EVD <- eigen(Re(gg))
                angles <- EVD$vectors[,1]
                if (angles[1]<0) angles <- -angles
                CC <- sqrt(0.5 + angles[1]/2)
                SS <- 0.5 * (complex(real=angles[2], imaginary=-angles[3]))/CC

                if(abs(SS)>eps){
                    encore <- TRUE
                    pair <- c(i, j)
                    G <- matrix(c(CC, -Conj(SS), SS, CC), nrow=2, byrow=TRUE)
                    V[,pair] <- V[,pair]%*%G
                    X[pair,] <- Conj(t(G)) %*% X[pair,]
                    X[,c(Ii,Ij)] <-  cbind(CC*X[,Ii]+SS*X[,Ij], -Conj(SS)*X[,Ii]+CC*X[,Ij])
                    }
                }
            }
        if (iter >= maxiter) stop("maxiter reached without convergence")
        }
    D<-array(X,c(p,p,k))
    if(type == "Matrix"){
        D <- aperm(D, c(1,3,2))
        D <- matrix(D, ncol=p)
        }
    RES <- list(V=V, D=D)
    RES
    }
