## sparse matrix version function for B'B
## default is "spam" if fed a listw of style "W" or "S", "naive" if a matrix
##
## usage and computational gains on large example:
#> system.time(b1s<-solveB(lambda=0.2, listw=lwUScw))
#   user  system elapsed
#   0.04    0.00    0.16
#> system.time(b1m<-solveB(lambda=0.2, listw=lwUScw, method="Matrix"))
#   user  system elapsed
#   0.14    0.00    0.55
#> system.time(b1n<-solveB(lambda=0.2, listw=lwUScw, method="naive"))
#   user  system elapsed
#  20.75    0.14   20.91

## needed here, but already in 'splm'??

listw2U_spam <- function(lw){
0.5 * (lw + t(lw))
}

similar.listw_spam <- function(listw)
{
    nbsym <- attr(listw$neighbours, "sym")
    if (is.null(nbsym))
        nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
    if (!nbsym)
        stop("Only symmetric nb can yield similar to symmetric weights")
    if (attr(listw$weights, "mode") == "general")
        if (!attr(listw$weights, "glistsym"))
            stop("General weights must be symmetric")
    n <- length(listw$neighbours)
    if (n < 1)
        stop("non-positive number of entities")
    sww <- as.spam.listw(listw)
    if (listw$style == "W") {
        sd <- attr(listw$weights, "comp")$d
        sd1 <- 1/(sqrt(sd))
        sdd <- diag.spam(sd, n, n)
        sdd1 <- diag.spam(sd1, n, n)
        sww1 <- sdd %*% sww
        res <- sdd1 %*% sww1 %*% sdd1
    }
    else if (listw$style == "S") {
        q <- attr(listw$weights, "comp")$q
        Q <- attr(listw$weights, "comp")$Q
        eff.n <- attr(listw$weights, "comp")$eff.n
        q1 <- 1/(sqrt(q))
        qq <- diag.spam(q, n, n)
        qq1 <- diag.spam(q1, n, n)
        ww0 <- (Q/eff.n) * sww
        ww1 <- qq %*% ww0
        sim0 <- qq1 %*% ww1 %*% qq1
        res <- (eff.n/Q) * sim0
    }
    else stop("Conversion not suitable for this weights style")
    res
}

xprodB <- function(lambda, listw, can.sim=TRUE) {

    ## check if listw or matrix;
    if("listw" %in% class(listw)) {
         ## case listw is 'listw'
         n <- length(listw$weights)
         if (listw$style %in% c("W", "S") & can.sim) {
            csrw <- listw2U_spam(similar.listw_spam(listw))
            similar <- TRUE
         } else {
             csrw <- as.spam.listw(listw)
         }
      } else {
         ## case listw is 'matrix'
         n <- ncol(listw)
         csrw <- as.spam(listw)
     }
     I <- diag.spam(1, n, n)
     B <-  I - lambda * csrw
     xprodB <- t(B) %*% B
     return(xprodB)
}

ldetB <- function(lambda, listw, can.sim=TRUE) {

    ## check if listw or matrix;
    if("listw" %in% class(listw)) {
        ## case listw is 'listw'
        if (listw$style %in% c("W", "S") & can.sim) {
            csrw <- listw2U_spam(similar.listw_spam(listw))
            similar <- TRUE
        }
        else csrw <- as.spam.listw(listw)
        n <- length(listw$weights)
        I <- diag.spam(1, n, n)
        B <- I - lambda * csrw
        J1 <- try(determinant(B, logarithm = TRUE)$modulus, silent = TRUE)
        if (class(J1) == "try-error") {
            ## fall back on standard methods
            ldetB <- log(det(as.matrix(B)))
            warning("Bad result in calculating log(det(B))")
        } else {
            ldetB <- J1
        }
     } else {
         ## case listw is 'matrix'
         n <- ncol(listw)
         ldetB <- log(det(diag(1, n) - lambda * listw))
         #csrw <- as.spam(listw)
     }

     return(ldetB)
}


