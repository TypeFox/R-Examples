context("Utility functions")

test_that("By", {
    b1 <- By(datasets::CO2,~Treatment+Type,colMeans,~conc)
    b2 <- By(datasets::CO2,c('Treatment','Type'),colMeans,'conc')
    expect_equivalent(b1,b2)
    ## require('data.table')
    ## t1 <- as.data.frame(data.table(datasets::CO2)[,mean(uptake),by=.(Treatment,Type,conc>500)])
    d0 <- transform(datasets::CO2,conc500=conc>500)
    t1 <- by(d0[,"uptake"],d0[,c("Treatment","Type","conc500")],mean)
    t2 <- By(datasets::CO2,~Treatment+Type+I(conc>500),colMeans,~uptake)
    expect_true(inherits(t2,"array"))
    expect_equivalent(sort(t2),sort(t1))
})


test_that("Expand", {
    dd <- Expand(iris, Sepal.Length=2:8, Species=c("virginica","setosa"))
    expect_identical(levels(iris$Species),levels(dd$Species))
    expect_true(nrow(dd)==14)
    
    d0 <- datasets::warpbreaks[,c("wool","tension")]
    T <- table(d0)
    d1 <- Expand(T)
    expect_identical(dim(d0),dim(d1))
    expect_identical(table(d1),T)

    expect_identical(expand.grid(1:2,1:2),Expand(1:2,1:2))
    expect_identical(expand.grid(a=1:2,b=1:2),Expand(a=1:2,b=1:2))
})

test_that("dsort", {
    data(hubble)
    expect_equivalent(order(dsort(hubble, ~sigma)$sigma),
                      seq_len(nrow(hubble)))
})

test_that("formulas", {
    f <- toformula(c('y1','y2'),'x'%++%1:5)
    ff <- getoutcome(f)
    expect_equivalent(trim(ff,all=TRUE),"c(y1,y2)")
    expect_true(length(attr(ff,'x'))==5)
})

test_that("trim", {
    expect_true(length(grep(" ",trim(" test ")))==0)    
    expect_true(length(gregexpr(" ",trim(" t e s t "))[[1]])==3)
    expect_true(length(grep(" ",trim(" t e s t ",all=TRUE)))==0)
})


test_that("Matrix operations:", {
    ## vec operator
    expect_equivalent(vec(diag(3)),c(1,0,0,0,1,0,0,0,1))
    expect_true(nrow(vec(diag(3),matrix=TRUE))==9)

    ## commutaion matrix
    A <- matrix(1:16 ,ncol=4)
    K <- commutation(A)
    expect_equivalent(K%*%as.vector(A),vec(t(A),matrix=TRUE))

    ## Block diagonal
    A <- diag(3)+1
    B <- blockdiag(A,A,A,pad=NA)
    expect_equivalent(dim(B),c(9,9))
    expect_true(sum(is.na(B))==81-27)
})


test_that("plotConf", {
    m <- lvm(y~x+g)
    distribution(m,~g) <- binomial.lvm()
    d <- sim(m,50)
    l <- lm(y~x+g,d)
    g1 <- plotConf(l,var2="g",plot=FALSE)
    g2 <- plotConf(l,var1=NULL,var2="g",plot=FALSE)
})


test_that("wrapvev", {
    expect_equivalent(wrapvec(5,2),c(3,4,5,1,2))
    expect_equivalent(wrapvec(seq(1:5),-1),c(5,1,2,3,4))
})


test_that("matrix functions", {
    A <- revdiag(1:3)
    expect_equivalent(A,matrix(c(0,0,1,0,2,0,3,0,0),3))
    expect_equivalent(1:3,revdiag(A))
    revdiag(A) <- 4
    expect_equivalent(rep(4,3),revdiag(A))
    diag(A) <- 0
    offdiag(A) <- 5
    expect_true(sum(offdiag(A))==6*5)

    
    A <- matrix(0,3,3)
    offdiag(A,type=3) <- 1:6
    B <- crossprod(A)
    
    expect_equivalent(solve(A),Inverse(A))
    expect_equivalent(det(B),attr(Inverse(B,chol=TRUE),"det"))
})


test_that("getmplus", {
})
test_that("getsas", {
})



test_that("All the rest", {
    expect_false(lava:::versioncheck(NULL))
    expect_true(lava:::versioncheck("lava",c(1,4,1)))

    op <- lava.options(debug=TRUE)
    expect_true(lava.options()$debug)
    lava.options(op)

    A <- diag(2); colnames(A) <- c("a","b")    
    expect_output(lava:::printmany(A,A,2,rownames=c("A","B"),bothrows=FALSE),"a b")
    expect_output(lava:::printmany(A,A[1,,drop=FALSE],2,rownames=c("A","B"),bothrows=FALSE),"a b")
    expect_output(printmany(A,A,2,rownames=c("A","B"),name1="no.1",name2="no.2",
                            bothrows=TRUE),"no.1")
    
})


