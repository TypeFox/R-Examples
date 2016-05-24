ReedfrogSizepred <- 
  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
             Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L))

rfs_dat <- c(list(nobs=nrow(ReedfrogSizepred),
                nexposed=rep(10,nrow(ReedfrogSizepred))),
                ReedfrogSizepred)
rfs_params <- list(c=0.45,d=13,g=1)
library(R2admb)
setup_admb()
m1 <- do_admb("ReedfrogSizepred0",
              data=rfs_dat,
              params=rfs_params,
              run.opts=run.control(checkparam="write",
                checkdata="write",clean=FALSE),
              verbose=TRUE)

proftime <- system.time(m1P <- do_admb("ReedfrogSizepred0",
              data=rfs_dat,
               params=rfs_params,
               profile=TRUE,
               profpars=c("c","d","g"),
               run.opts=run.control(checkparam="write",
                 checkdata="write"),
               admb_errors="warn",  ## because of Hessian problem
               verbose=TRUE))
mctime <- system.time(m1MC <- do_admb("ReedfrogSizepred0",
                                     data=rfs_dat,
                                      params=rfs_params,
                                      run.opts=run.control(checkparam="write",
                                        checkdata="write",clean=FALSE),
                                      mcmc=TRUE,
                                      mcmc.opts=mcmc.control(mcmcpars=c("c","d","g")),
                                      verbose=TRUE))
m1MC2 <- do_admb("ReedfrogSizepred0",
                 data=rfs_dat,
                 params=rfs_params,
                 run.opts=run.control(checkparam="write",
                   checkdata="write"),
                 mcmc=TRUE,
                 mcmc.opts=mcmc.control(mcmcpars=c("c","d","g"),
                   mcmc=5000),
                verbose=TRUE)
## ???  something a bit funny going on here but can't quite figure it out
## run_admb("ReedfrogSizepred0_gen",verbose=TRUE)
## run_admb("ReedfrogSizepred0_gen",mcmc=TRUE,verbose=TRUE)
## run_admb("ReedfrogSizepred0_gen",mcmc=TRUE,
##          mcmc.opts=mcmc.control(mcmc=5000),verbose=TRUE)

unlink("reedfrogsizepred0")
save("m1","m1MC","m1P","m1MC2","mctime","proftime",file="Reedfrog_runs.RData")

