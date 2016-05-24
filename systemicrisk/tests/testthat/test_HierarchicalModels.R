context("Test of Hierarchical Models")

test_that("sample_HierachicalModel test presence of nondefault enforced 0",{
    skip_on_cran()

    set.seed(12930)
    n <- 4
    replicate(10,{
        p <- matrix(nrow=n,ncol=n,sample(seq(0,1,length.out=n*n)))
        mp <- Model.p.constant(n=n,p=p)
        m <- Model.Indep.p.lambda(model.p=mp,
                                  model.lambda=Model.lambda.constant(n=n,lambda=1/2))
        L <- genL(m)$L
        res <- sample_HierarchicalModel(l=rowSums(L),a=colSums(L),model=m,nsamples=1,thin=100,burnin=2,silent=TRUE)
        expect_true(all(ifelse(p==0,res$L[[1]]==0,TRUE)))
    })
})
