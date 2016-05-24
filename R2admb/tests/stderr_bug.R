## NOT RUN: won't work on CRAN
## FIXME: set up environment variables to allow this to run
## if not on CRAN
library("R2admb")
tadpoledat <-
    data.frame(TBL = rep(c(9,12,21,25,37),each=3),
               Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L),
               nexposed=rep(10,15))

if (FALSE) {
    setup_admb()
    file.copy(system.file("tplfiles",
                          "ReedfrogSizepred0.tpl",package="R2admb"),
              "tadpole.tpl")
     m1 <- do_admb("tadpole",
                   data=c(list(nobs=15),tadpoledat),
                   params=list(c=0.45,d=13,g=1),
                   bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
                   run.opts=run.control(checkparam="write",
                   checkdata="write",clean=FALSE))
    unlink("tadpole.tpl")
    summary(m1)
    read.table("tadpole_gen.std",header=TRUE)[,4]
    sqrt(diag(read_admb("tadpole_gen")$vcov))
    sqrt(diag(read_pars("tadpole_gen")$vcov))
    ## compare with mle2 results:
    ## library("bbmle")
    ## m2 <- mle2(Kill~dbinom(prob=c*(TBL/d*exp(1-TBL/d))^g,size=nexposed),
    ##      data=tadpoledat,
    ##      start=as.list(coef(m1)))
    m2_vcov <- structure(c(0.0158042534788236,
                           0.0578065701530253, 0.504366729952002, 
                           0.0578065701530253, 0.657825248678555,
                           2.24662255486992, 0.504366729952002, 
                           2.24662255486992, 36.3935910094069),
                         .Dim = c(3L, 3L), .Dimnames = list(
                                           c("c", "d", "g"),
                                           c("c", "d", "g")))
    all.equal(vcov(m1),m2_vcov,tol=2e-4)
}
