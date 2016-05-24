
FIM <- evaluate.fim(poped.db) 
dmf <- det(FIM)


blockfinal(fn="",fmf=FIM,
           dmf=dmf,
           groupsize=poped.db$design$groupsize,
           ni=poped.db$design$ni,
           xt=poped.db$design$xt,
           x=poped.db$design$x,a=poped.db$design$a,
           model_switch=poped.db$design$model_switch,
           poped.db$parameters$param.pt.val$bpop,
           poped.db$parameters$param.pt.val$d,
           poped.db$parameters$docc,
           poped.db$parameters$param.pt.val$sigma,
           poped.db,
           opt_xt=TRUE,
           fmf_init=FIM,
           dmf_init=dmf,
           param_cvs_init=rbind(get_rse(FIM,poped.db,use_percent=FALSE)))


