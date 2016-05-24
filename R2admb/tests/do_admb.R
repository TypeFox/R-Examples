library(R2admb)

s <- try(setup_admb())
if (!identical(aa <- admb_version(),NA)) {
    ## run only if we can
    file.copy(system.file("tplfiles","ReedfrogSizepred0.tpl",package="R2admb"),"tadpole.tpl")
    tadpoledat <-  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
                              Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L),
                              nexposed=rep(10,15))
    par1 <- list(c=0.45,d=13,g=1)
    ## getvals(par1,"value",valsOK=TRUE)
    m1 <- do_admb("tadpole",
                  data=c(list(nobs=15),tadpoledat),
                  params=par1,
                  bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
                  run.opts=run.control(checkparam="write",
                  checkdata="write",clean="all"))

    m1P <- do_admb("tadpole",
                   data=c(list(nobs=15),tadpoledat),
                   params=par1,
                   bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
                   phase=list(c=1,d=2,g=-1),
                   run.opts=run.control(checkparam="write",
                   checkdata="write",clean="all"))

    par2 <- list(c=list(0.45,bounds=c(0,1)),
                 d=list(13,bounds=c(0,50)),
                 g=list(1,bounds=c(-1,25)))

    m2 <- do_admb("tadpole",
                  data=c(list(nobs=15),tadpoledat),
                  params=par2,
                  run.opts=run.control(checkparam="write",
                  checkdata="write",clean="all"))

    par2P <- list(c=list(0.45,bounds=c(0,1),phase=1),
                 d=list(13,bounds=c(0,50),phase=2),
                 g=list(1,bounds=c(-1,25),phase=-1))

    m2P <- do_admb("tadpole",
                  data=c(list(nobs=15),tadpoledat),
                  params=par2P,
                  run.opts=run.control(checkparam="write",
                  checkdata="write",clean="all"))

    stopifnot(all.equal(m1,m2))
    stopifnot(all.equal(m1P,m2P))
    unlink("tadpole.tpl")

    file.copy(system.file("tplfiles","toy1.tpl",package="R2admb"),"toy1.tpl")
    ## random effects example??
    cbpp <-  structure(list(herd = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 
         3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 
         7L, 7L, 7L, 7L, 8L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 11L, 
         11L, 11L, 11L, 12L, 12L, 12L, 12L, 13L, 13L, 13L, 13L, 14L, 14L, 
         14L, 14L, 15L, 15L, 15L, 15L),
     .Label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"),
                            class = "factor"), 
    incidence = c(2, 3, 4, 0, 3, 1, 1, 8, 2, 0, 2, 2, 0, 2, 0, 
    5, 0, 0, 1, 3, 0, 0, 1, 8, 1, 3, 0, 12, 2, 0, 0, 0, 1, 1, 
    0, 2, 0, 5, 3, 1, 2, 1, 0, 0, 1, 2, 0, 0, 11, 0, 0, 0, 1, 
    1, 1, 0), size = c(14, 12, 9, 5, 22, 18, 21, 22, 16, 16, 
    20, 10, 10, 9, 6, 18, 25, 24, 4, 17, 17, 18, 20, 16, 10, 
    9, 5, 34, 9, 6, 8, 6, 22, 22, 18, 22, 25, 27, 22, 22, 10, 
    8, 6, 5, 21, 24, 19, 23, 19, 2, 3, 2, 19, 15, 15, 15), period = structure(c(1L, 
    2L, 3L, 4L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 
    2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 1L, 2L, 3L, 
    4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 
    3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L), .Label = c("1", 
    "2", "3", "4"), class = "factor")), .Names = c("herd", "incidence", 
"size", "period"), row.names = c(NA, 56L), class = "data.frame")

    X <- model.matrix(~period,data=cbpp)
    Zherd <- model.matrix(~herd-1,data=cbpp)

    tmpdat <- list(X=X,Zherd=Zherd,
                 incidence=cbpp$incidence,size=cbpp$size,
                 nobs=nrow(cbpp))

    d1 <- do_admb("toy1",
                  data=tmpdat,
              params=list(beta=rep(0,ncol(X)),sigma_herd=0.1),
              bounds=list(sigma_herd=c(0.0001,20)),
              re=list(u_herd=ncol(Zherd)),
              run.opts=run.control(checkdata="write",checkparam="write"),
              mcmc=TRUE,
              mcmc.opts=mcmc.control(mcmc=20,mcmcpars=c("beta","sigma_herd")))

        d2 <- do_admb("toy1",
                  data=tmpdat,
              params=list(beta=rep(0,ncol(X)),
                  sigma_herd=list(0.1,bounds=c(0.0001,20))),
              re=list(u_herd=ncol(Zherd)),
              run.opts=run.control(checkdata="write",checkparam="write"))

    summary(d2)
    unlink("toy1.tpl")
}
