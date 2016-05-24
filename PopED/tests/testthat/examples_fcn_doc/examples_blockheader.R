
FIM <- evaluate.fim(poped.db) 
dmf <- det(FIM)

blockheader(poped.db,name="")

blockheader(name="",iter=1,poped.db)


blockheader(name='',
              iter=1,
              poped.db,
              e_flag=FALSE,
              opt_xt=TRUE,
              opt_a=TRUE,opt_x=poped.db$settings$optsw[4],
              opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
              fmf=FIM,dmf=dmf,
              bpop=poped.db$parameters$param.pt.val$bpop,
              d=poped.db$parameters$param.pt.val$d,
              docc=poped.db$parameters$docc,sigma=poped.db$parameters$param.pt.val$sigma)



blockheader(name='',
              iter=1,
              poped.db,
              e_flag=TRUE,
              opt_xt=TRUE,
              opt_a=TRUE,opt_x=poped.db$settings$optsw[4],
              opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
              fmf=FIM,dmf=dmf,
              bpop=poped.db$parameters$param.pt.val$bpop,
              d=poped.db$parameters$param.pt.val$d,
              docc=poped.db$parameters$docc,sigma=poped.db$parameters$param.pt.val$sigma)
  
  
poped.db.1 <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                  fg_file="sfg",
                                  fError_file="feps.add.prop",
                                  bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                  notfixed_bpop=c(1,1,1,0),
                                  d=c(CL=0.07, V=0.02, KA=0.6), 
                                  sigma=c(0.01,0.25),
                                  groupsize=32,
                                  xt=rbind(c( 0.5,1,2,6,24,36,72,120),
                                           c( 0.5,1.1,2,6,24,36,72,120)),
                                  minxt=rbind(c(0,1,1.5,3,20,30,70,118),
                                              c(0.1,1.1,1.6,3.1,20.1,30.1,70.1,118.1)),
                                  maxxt=c(12,13,14,15,26,44,78,120),
                                  a=70,
                                  mina=0,
                                  maxa=100)


blockheader(poped.db.1,name="",trflag=2,opt_xt=TRUE)


