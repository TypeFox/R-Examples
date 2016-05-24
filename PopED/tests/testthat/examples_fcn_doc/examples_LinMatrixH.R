
#for the FO approximation
ind=1
LinMatrixH(model_switch=t(poped.db$design$model_switch[ind,,drop=FALSE]),
          xt_ind=t(poped.db$design$xt[ind,,drop=FALSE]),
          x=zeros(0,1),
          a=t(poped.db$design$a[ind,,drop=FALSE]),
          bpop=poped.db$parameters$bpop[,2,drop=FALSE],
          b_ind=zeros(poped.db$parameters$NumRanEff,1),
          bocc_ind=zeros(poped.db$parameters$NumDocc,1),
          poped.db)[["y"]]


