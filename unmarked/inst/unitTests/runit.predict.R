
test.occu <- function() {

    if(!require(raster))
        stop("raster package required")
    set.seed(55)
    R <- 20
    J <- 4
    x1 <- rnorm(R)
    x2 <- factor(c(rep('A', R/2), rep('B', R/2)))
    x3 <- matrix(rnorm(R*J), R, J)
    z <- rbinom(R, 1, 0.5)
    y <- matrix(rbinom(R*J, 1, z*0.6), R, J)
    x1[1] <- NA
    x3[2,1] <- NA
    x3[3,] <- NA
    umf1 <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x1=x1, x2=x2),
                              obsCovs=list(x3=x3))
    fm1 <- occu(~x3 ~x1+x2, umf1)
    E1.1 <- predict(fm1, type="state")
    E1.2 <- predict(fm1, type="det")

    nd1.1 <- data.frame(x1=0, x2=factor('A', levels=c('A','B')))
    nd1.2 <- data.frame(x3=0)
    E1.3 <- predict(fm1, type="state", newdata=nd1.1, appendData=TRUE)
    E1.4 <- predict(fm1, type="det", newdata=nd1.2)

    r1 <- raster(matrix(rnorm(100), 10))
    checkException(predict(fm1, type="state", newdata=r1))
    s1 <- stack(r1)
    checkException(predict(fm1, type="state", newdata=s1))
    names(s1) <- c("x3")
    E1.5 <- predict(fm1, type="det", newdata=s1)
    E1.5 <- predict(fm1, type="det", newdata=s1, appendData=TRUE)

    E1.6 <- predict(fm1, type="state", level=0.9)
    checkEquals(as.numeric(E1.6[1,3:4]), c(0.01881844, 0.8538048))

}



test.pcount <- function() {

    set.seed(55)
    R <- 20
    J <- 4
    N <- rpois(R, 2)
    y <- matrix(rbinom(R*J, N, 0.7), R, J)
    umf1 <- unmarkedFramePCount(y=y)

    fm1 <- pcount(~1 ~1, umf1, K=40)
    E1.1 <- predict(fm1, type="state")
    E1.2 <- predict(fm1, type="det")

    fm2 <- pcount(~1 ~1, umf1, K=40, mixture="NB")
    E2.1 <- predict(fm2, type="state")
    checkException(predict(fm2, type="alpha"))

    fm3 <- pcount(~1 ~1, umf1, K=40, mixture="ZIP")
    E3.1 <- predict(fm3, type="state")
    checkException(predict(fm3, type="psi"))
    checkEquals(E3.1[1,1], 1.818512, tol=1e-6)

}





test.pcountOpen <- function() {

    set.seed(55)
    R <- 10
    J <- 4
    T <- 3
    N <- matrix(NA, R, T)
    N[,1] <- rpois(R, 4)
    N[,2] <- rbinom(R, N[,1], 0.8) + rpois(R, 1)
    N[,3] <- rbinom(R, N[,2], 0.8) + rpois(R, 1)
    y1 <- matrix(rbinom(R*J, N[,1], 0.7), R, J)
    y2 <- matrix(rbinom(R*J, N[,2], 0.7), R, J)
    y3 <- matrix(rbinom(R*J, N[,3], 0.7), R, J)
    umf1 <- unmarkedFramePCO(y=cbind(y1,y2,y3), numPrimary=T)

#    fm1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=30)
#    E1.1 <- predict(fm1, type="lambda")
#    E1.2 <- predict(fm1, type="det")

    fm2 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=40, mixture="NB")
    checkException(predict(fm2, type="alpha"))

    fm3 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=40, mixture="ZIP")
    E3.1 <- predict(fm3, type="lambda")
    checkException(predict(fm3, type="psi"))
    checkEquals(E3.1[1,1], 2.472029, tol=1e-6)

}



