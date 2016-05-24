context("getcdfupdates")

testmodel <- function(model,data,testpoints=seq(-5,5,length.out=201)){
    P <- model$Pofdata(data)
    xi <- model$xiofP(P)
    Fcadlag <- model$getcdfupdates(P,xi,cadlag=TRUE)
    Fcaglad <- model$getcdfupdates(P,xi,cadlag=FALSE)
    expect_true(all(Fcadlag(testpoints)>=0))
    expect_true(all(Fcadlag(testpoints)<=1))
    expect_true(all(Fcaglad(testpoints)>=0))
    expect_true(all(Fcaglad(testpoints)<=1))
    expect_true(all(Fcadlag(testpoints)>=Fcaglad(testpoints)-1e-15))

}


test_that("SPCModelNormal  equality in getcdfupdate",{
    for (Delta in c(-1,0,2,5)){
        testmodel(SPCModelNormal(Delta=Delta),data=rnorm(100))
        testmodel(SPCModelNormal(Delta=Delta),data=seq(-5,5,length.out=201))
        testmodel(SPCModelNormal(Delta=Delta),data=seq(-5,5,length.out=201))
    }
})


test_that("SPCModelNonpar  equality in getcdfupdate",{
    testmodel(SPCModelNonpar(update=function(xi,data) data, xiofP=function(P)list()),
              data=rnorm(100))
    testmodel(SPCModelNonpar(update=function(xi,data) data, xiofP=function(P)list()),
              data=seq(-5,5,length.out=201))
    testmodel(SPCModelNonpar(update=function(xi,data) data, xiofP=function(P)list()),
              data=sample(seq(-5,5,length.out=201)))
})
