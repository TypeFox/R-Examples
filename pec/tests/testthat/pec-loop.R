test_that("boot632plus",{
    set.seed(130971)
    dat <- SimSurv(100)
    dat$X1 <- as.factor(dat$X1)
    Models <- list("Cox.X1"=coxph(Surv(time,status)~X1,data=dat,y=TRUE),
                   "Cox.X2"=coxph(Surv(time,status)~X2,data=dat,y=TRUE),
                   "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=dat,y=TRUE))
    set.seed(17100)
    PredError.632plus <- pec(object=Models,
                             formula=Surv(time,status)~X1+X2,
                             data=dat,
                             exact=TRUE,
                             cens.model="cox",
                             splitMethod="Boot632plus",
                             B=3,
                             verbose=TRUE)
    saved <- c(Reference=8.09,Cox.X1=7.86,Cox.X2=7.23,Cox.X1.X2=7.07)
    expect_equal(round(100*ibs(PredError.632plus,times=3)[,4],2),saved)
})
