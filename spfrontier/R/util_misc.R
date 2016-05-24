# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Miscellaneous utils
#


paramsToVector <- function (p, olsen = T) {
    vec <- c()
    if(olsen) {
        vec <- c(p$gamma)
        vecNames <- paste("Gamma", seq(length(p$gamma)), sep = "")
    }else{
        vec <- c(p$beta)
        vecNames <- paste("Beta", seq(length(p$beta)), sep = "")
    }
    if (!is.null(p$rhoY)){
        vec <- c(vec, p$rhoY)
        vecNames <- c(vecNames, "rhoY")
    }
    if(olsen) {
        vec <- c(vec, p$nu, p$lambda)
        vecNames <- c(vecNames, "Nu", "Lambda")
    }else{
        vec <- c(vec, p$sigmaV, p$sigmaU)
        vecNames <- c(vecNames, "sigmaV", "sigmaU")
    }
    if (!is.null(p$rhoV)){
        vec <- c(vec, p$rhoV)
        vecNames <- c(vecNames, "rhoV")
    }
    if (!is.null(p$rhoU)){
        vec <- c(vec, p$rhoU)
        vecNames <- c(vecNames, "rhoU")
    }
    if (!is.null(p$mu)){
        vec <- c(vec, p$mu)
        vecNames <- c(vecNames, "mu")
    }
    names(vec) <- vecNames
    return(vec)
}

paramsFromVector <- function (parameters, k, isSpY, isSpV, isSpU, isTN, olsen = T) {
    ret <- list()
    if (olsen) {
        ret$gamma <- parameters[1:k]
    }else{
        ret$beta <- parameters[1:k]
    }
    p <- k+1
    rhoY <- NULL
    if (isSpY){
        ret$rhoY <- parameters[p]
        p <- p+1
    }
    if (olsen) {
        ret$nu <- parameters[p]
        ret$lambda <- parameters[p+1]
    }else{
        ret$sigmaV <- parameters[p]
        ret$sigmaU <- parameters[p+1]
    }
    p <- p + 2
    if (isSpV){
        ret$rhoV <- parameters[p]
        p <- p + 1
    }
    if (isSpU){
        ret$rhoU <- parameters[p]
        p <- p + 1
    }
    if(isTN){
        ret$mu <- parameters[p]
    }
    return(ret)
}

olsenReparamBack = function(par){
    sigma <- 1/par$nu
    par$beta <- as.vector(par$gamma*sigma)
    par$sigmaV <- as.numeric(sigma/sqrt(1+par$lambda^2))
    par$sigmaU <- as.numeric(par$lambda*par$sigmaV)
    
    par$gamma <- NULL
    par$nu <- NULL
    par$lambda <- NULL
    
    X <- envirGet("X")
    names(par$beta) <- colnames(X)
    names(par$sigmaV) <- c("sigmaV")
    names(par$sigmaU) <- c("sigmaU")
    return(par)
}

olsenReparam = function(par){
    sigma <- sqrt(par$sigmaV^2+par$sigmaU^2)
    par$lambda <- par$sigmaU/par$sigmaV
    par$nu <- 1/sigma
    par$gamma <- par$beta/sigma
    
    par$beta <- NULL
    par$sigmaU <- NULL
    par$sigmaV <- NULL
    
    names(par$gamma) = paste("Gamma", seq(length(par$gamma)), sep = "")
    names(par$nu) <- c("nu")
    names(par$lambda) <- c("lambda")
    return(par)
}

# olsenGradient<- function(par){
#     p <- paramsToVector(par, olsen=F)
#     G <- diag(length(p))
#     sigma <- sqrt(par$sigmaV^2+par$sigmaU^2)
#     b <- length(par$beta)
#     G[1:b,1:b] <- diag(b)/sigma
#     j <- b+1
#     if (!is.null(par$rhoY)){
#         j <- j+1
#     }
#     G[1:b,j:j] <- -par$beta*par$sigmaV/sigma^3
#     G[j,j] <- -par$sigmaV/sigma^3
#     G[(j+1),j] <- -par$sigmaU/par$sigmaV^2
#     j <- j+1
#     G[1:b,j:j] <- -par$beta*par$sigmaU/sigma^3
#     G[(j-1),j] <- -par$sigmaU/sigma^3
#     G[j,j] <- 1/par$sigmaV
#     return(G)
# }

olsenGradient<- function(par){
    p <- paramsToVector(par, olsen=T)
    
    G <- diag(length(p))
    b <- length(par$gamma)
    G[1:b,1:b] <- diag(b)/par$nu
    j <- b+1
    if (!is.null(par$rhoY)){
        j <- j+1
    }
    a <- sqrt(1+par$lambda^2)
    G[1:b,j:j] <- -par$gamma/par$nu^2
    G[j,j] <- -1/(a*par$nu^2)
    G[(j+1),j] <- -par$lambda/(a*par$nu^2)
    j <- j+1
    G[(j-1),j] <- -par$lambda/(a^3*par$nu)
    G[j,j] <- 1/(a^3*par$nu)
    return(G)
}



#' @title Standard spatial contiguity matrixes
#'
#' @description
#' \code{genW} generates an spatial contiguity matrix (rook or queen)
#' 
#' @details
#' To generate spatial interaction between \code{n} objects the function arranges them on a chess board. 
#' A number of columns is calculated as a square root of \code{n}, rounded to the top. The last row contains empty cells, if n is not quadratic 
#' 
#' 
#' @param n a number of objects with spatial interaction to be arranged.See 'Details' for objects arranging principle
#' @param type an optional type of spatial interaction. Currently 'rook' and ''queen' values are supported, to produce Rook and Queen Contiguity matrix. See references for more info. 
#' By default set to rook.
#' @param seed an optional random number generator seed for random matrices 
#' 
#' @keywords spatial stochastic frontier
#' @references 
#' Anselin, L. (1988). Spatial Econometrics: Methods and Models. Kluwer Academic Publishers, Dordrecht, The Netherlands.
#' 
#' @rdname util-misc
#' @examples
#' # Completely filled 10x10 rook contiguity matrix
#' rookW <- genW(100)
#' rookW
#' 
#' # Partly filled 10x10 rook contiguity matrix
#' rookW <- genW(90)
#' rookW
#' 
#' # Completely filled 10x10 queen contiguity matrix
#' queenW <- genW(100, type="queen")
#' queenW

genW <- function(n, type="rook", seed=NULL){
    W <- matrix(nrow=n, ncol=n, rep(0,n*n))
    if (type=="rook" || type=="queen"){
        k <- ceiling(sqrt(n))
        for (i in 1:n){
            h <- (i-1) %/% k + 1
            v <- (i-1) %% k + 1
                for (dv in c(-1,0,1)){
                    for (dh in c(-1,0,1)){
                        diff <- abs(dv) + abs(dh)
                        v1 <- v + dv
                        h1 <- h + dh
                        neib <- (diff>0) & (v1 >0) & (h1 >0) & (v1 <=k) & (h1 <=k) & (type=="queen" || diff==1)
                        if(neib){
                            j <- v1 + (h1-1)*k
                            if (j<=n) W[i,j] = 1
                        }
                    }
                }
        }
    }else if (type=="inverse distance"){
        if (!is.null(seed)) set.seed(seed);
        coords <- matrix(nrow=n, ncol=2, rnorm(n*2, mean=0, sd=1))
        #plot(coords)
        W <- constructW(coords, seq(1, n))
        W[which(W<quantile(W,probs=0.5))] <- 0
    }
    return(W)
}



#' @title Spatial contiguity matrix standartisation
#'
#' @description
#' \code{rowStdrt} standartizes spatial contiguity matrix by rows
#' 
#' @details
#' The function divides every element in an argument matrix by the sum of elements in its row. Some spatial estimation requires this standartisation (generally - for faster calculations)
#' 
#' 
#' @param W a spatial contiguity matrix to be standatised
#' 
#' @rdname util-misc
#' @examples
#' 
#' # Completely filled 10x10 queen contiguity matrix
#' queenW <- genW(100, type="queen")
#' queenW
#' 
#' # Standartisation
#' stQueenW <- rowStdrt(queenW)
#' stQueenW
#' 

rowStdrt = function(W){
    for (j in 1:nrow(W)){
        W[j,] <- W[j,]/sum(W[j,])
    }
    return(W)
}


#' @title Spatial contiguity matrix construction
#'
#' @description
#' \code{constructW} contructs a spatial contiguity matrix using object longitude and latitude coordinates
#' 
#' @details
#' The function contructs a spatial contiguity matrix 
#' using object longitude and latitude coordinates. Eucledean distance is currently used.
#' 
#' 
#' @param coords a matrix of two columns, where every row is a longitude-latitude pair of object coordinates
#' @param labels a vector of object lables to mark rows and columns of the resulting contiguity matrix
#' 
#' @rdname util-misc
#' @examples
#' 
#' data(airports)
#' 
#' W <- constructW(cbind(airports$lon, airports$lat),airports$ICAO_code)
#' 

constructW <- function(coords, labels){
    W <- 1/as.matrix(dist(coords))
    colnames(W) <- labels
    rownames(W) <- labels
    W[which(W==Inf)] <- 0
    return(W)
}




#
# Returns significance code
#
pvalMark <- function(pval){
    m <- ""
    if (is.na(pval)){
        m <- "not defined"
    } else if (pval < 0.001){
        m <- "***"
    } else if (pval < 0.01){
        m <- "**"
    } else if (pval < 0.05){
        m <- "*"
    } else if (pval < 0.1){
        m <- "."
    }
    return(m)
}


