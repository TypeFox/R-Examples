context("Test Automatic Choice of Thinning")



test_that("Automatic Thinning 10 by 10",{
    skip_on_cran()

    n <- 10
    set.seed(12930)
    for (j in 1:3){
        m <- Model.Indep.p.lambda(Model.p.constant(n,p=0.8),
                                  Model.lambda.constant(n, lambda=1))
        x <- genL(m)
        nsamples<-1e3
        for (rESStarget in c(0.1,0.5,0.9)){
            thin <- choosethin(l=rowSums(x$L),a=colSums(x$L),model=m,
                               relESStarget=rESStarget,silent=TRUE)
            res <- sample_HierarchicalModel(l=rowSums(x$L),a=colSums(x$L),model=m,thin=thin,nsamples=nsamples,silent=TRUE)
            allL <- sapply(res$L,c)
            ress <- coda::effectiveSize(t(allL))
            ress <- ress[!as.logical(c(diag(rep(1,n))))]/nsamples
            expect_true(quantile(ress,0.05)>rESStarget*0.9)
        }
    }
})


test_that("Automatic Thinning BetaGammaPrior",{
    skip_on_cran()

    n <- 5
    set.seed(12930)
    for (j in 1:3){
        m <- Model.Indep.p.lambda(Model.p.BetaPrior(n,shape1=100,shape2=10),
                           Model.lambda.GammaPrior(n,scale=1e-1))

        x <- genL(m)
        nsamples<-1e3
        for (rESStarget in c(0.1,0.5,0.8)){
            thin <- choosethin(l=rowSums(x$L),a=colSums(x$L),model=m,
                               relESStarget=rESStarget,silent=TRUE)
            res <- sample_HierarchicalModel(l=rowSums(x$L),a=colSums(x$L),model=m,thin=thin,nsamples=nsamples,silent=TRUE)
            allL <- sapply(res$L,c)
            ress <- coda::effectiveSize(t(allL))
            ress <- ress[!as.logical(c(diag(rep(1,n))))]/nsamples
            expect_true(quantile(ress,0.05)>rESStarget*0.9)
        }
    }
})
