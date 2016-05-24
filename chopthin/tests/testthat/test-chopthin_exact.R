set.seed(1237123)

context("test-choptin_exact")

checkbasics <- function(res,N){
    expect_true(length(res$indices)==N);
    expect_true(length(res$weights)==N);
    expect_true(all(res$indices>0));
    expect_true(all(res$weights>=0));
    expect_true(abs(sum(res$weights)-N)<1e-10); #weights normalised to N
}

test_that("Individual small scale tests",{
    expect_error(chopthin(rep(0,5),5));
    checkbasics(chopthin(1,N=1),1)

    for (i in c(100,101,105,150,197,198,199,200)){
        res <- chopthin(rep(1,100),N=i);
        checkbasics(res,i)
     }
})


test_that("Eta=Inf tests",{
    expect_error(chopthin(c(1,2,3),N=4,eta=Inf));
    set.seed(123910)
    checkbasics(chopthin(c(rexp(50),rep(0,50)),N=100,eta=Inf),100)
    checkbasics(chopthin(c(1,0),N=2,eta=Inf),2)
    checkbasics(chopthin(rexp(100),N=50,eta=Inf),50)
    checkbasics(chopthin(rexp(100),N=100,eta=Inf),100)
    for (N in c(100,99,50,35,2,1)){
        w <- c(rexp(50),rep(0,50))
        checkbasics(chopthin(w,N=N,eta=Inf),N)
        expect_true(abs(sum(chopthin(w,N=N,eta=Inf,normalise=FALSE)$weights)-sum(w))<1e-10)
    }
})


test_that("samples",{
    set.seed(1230);
    w <- rexp(10);
    res <- replicate(1e4,chopthin(w,N=1)$indices)
    expect_true(chisq.test(tabulate(res,nbins=10),p=w/sum(w))$p.value>1e-1)

    skip_on_cran()
    for (i in 1:5){
        w <- rexp(20);
        res <- replicate(2e5,chopthin(w,N=1)$indices)
        tab <- tabulate(res,nbins=20)
        expect_true(chisq.test(tab,p=w/sum(w))$p.value>1e-3)
    }
})

test_that("Preservation of Weights",{
    set.seed(1293013)
    w <- rexp(10);
    for (i in 1:100){
        res <- chopthin(w,N=i,normalise=FALSE);
        expect_equal(sum(res$weights),sum(w))
    }
})

test_that("Preservation of Weights, eta=4, special case",{
    set.seed(1239019)
    w <- c(0.5,1.0,2.0,3.0,4.0)
    for (i in 1:10){
        res <- chopthin(w,N=6,eta=4,normalise=FALSE)
        expect_equal(sum(res$weights),sum(w))
        expect_true(max(res$weights)/min(res$weights)<=4)
    }
})


test_that("Ratio observed - special case, eta=4",{
    set.seed(1239019)
    w <- c(0.55,1.0,2.5,3.9,2.5)
    expect_true(max(replicate(1e3, {res <- chopthin(w,N=6,eta=4,normalise=FALSE); max(res$weights)/min(res$weights)}))<=4)
})

test_that("Particle to be chopped dies",{
    set.seed(121231)
    w <- c(0.5,1.0,2.0001,2.9999)
    expect_true(all(replicate(1e3, is.element(3,chopthin(w,N=4,eta=4,normalise=FALSE)$indices))))
})


test_that("Ratio not necessarily respected",{
    set.seed(129301)
    w <- c(0.5,1.0,4.4,2.6)
    expect_true(max(replicate(1e3,{res <- chopthin(w, N = 5, eta = 4, normalise = FALSE)$weights;max(res)/min(res)}))<4.0)
})
