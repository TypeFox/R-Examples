"lmomgpaRC" <-
function(para) {
    z <- list(lambdas = NULL,
              ratios  = NULL,
              source = "lmomgpaRC",
              message = "These are B-type L-moments",
              zeta = NULL
             )

    if(! are.pargpa.valid(para)) return()

    attributes(para$para) <- NULL

    XI <- para$para[1]
    A  <- para$para[2]
    K  <- para$para[3]


    zeta <- NULL
    if(is.na(para$zeta)) { zeta <- 1 } else { zeta <- para$zeta }

    mr  <- function(r,z,k) { (1 - (1-z)^(r+k))/(r+k) }
    mr1 <- mr(1,zeta,K)
    mr2 <- mr(2,zeta,K)
    mr3 <- mr(3,zeta,K)
    mr4 <- mr(4,zeta,K)
    mr5 <- mr(5,zeta,K)

    L <- vector(mode = "numeric", length=5)
    R <- vector(mode = "numeric", length=5)
    R[1] <- NA

    L[1] <- XI+A*mr1    # r = 1
    L[2] <- A*(mr1 -    mr2) # r = 2
    L[3] <- A*(mr1 -  3*mr2  +  2*mr3) # r = 3
    L[4] <- A*(mr1 -  6*mr2  + 10*mr3 -  5*mr4) # r = 4
    L[5] <- A*(mr1 - 10*mr2  + 30*mr3 - 35*mr4 +  14*mr5) # r = 5

    R[3] <- L[3]/L[2]
    R[4] <- L[4]/L[2]
    R[5] <- L[5]/L[2]
    R[2] <- L[2]/L[1]

    z$lambdas <- L
    z$ratios  <- R
    z$zeta    <- zeta
    return(z)
}


# THE FOLLOWING CODE CHUCK USED TO DETERMINE COEFFICIENTS ON
# Lambda COMPUTATIONS
#library(lmomco)

#zeta <- 1
#mr <- function(r,z,k) { (1 - (1-z)^(r+k))/(r+k) }
#K <- .1
#alpha <- 7
#xi <- 15
#    mr1 <- mr(1,zeta,K)
#    mr2 <- mr(2,zeta,K)
#    mr3 <- mr(3,zeta,K)
#    mr4 <- mr(4,zeta,K)
#    mr5 <- mr(5,zeta,K)

#para <- vec2par(c(xi,alpha,K),type='gpa')
#lmr <- lmomgpa(para)


#A <- seq(1,60,by=1)
#B <- seq(1,60,by=1)
#C <- seq(1,60,by=1)
#D <- seq(1,60,by=1)

#L5 <- lmr$L5

#mindiff <- Inf
#while(1) {
#  a <- A[as.integer(runif(1,min=1,max=length(A)))]
#  b <- A[as.integer(runif(1,min=1,max=length(B)))]
#  c <- A[as.integer(runif(1,min=1,max=length(C)))]
#  d <- A[as.integer(runif(1,min=1,max=length(D)))]
#  L5.hat <- alpha*(mr1 - a*mr2 + b*mr3 - c*mr4 + d*mr5)
#  diff <- abs(L5.hat - L5)
#  if(diff <= mindiff) {
#     mindiff <- diff
#     cat(c(diff,a,b,c,d,"\n"))
#  }
#}
