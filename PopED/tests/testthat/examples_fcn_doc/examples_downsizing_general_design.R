
dgd <- downsizing_general_design(poped.db)

output = mftot(dgd$model_switch,poped.db$design$groupsize,
               dgd$ni,dgd$xt,dgd$x,dgd$a,
               poped.db$parameters$param.pt.val$bpop,
               poped.db$parameters$param.pt.val$d,
               poped.db$parameters$param.pt.val$sigma,
               poped.db$parameters$param.pt.val$docc,
               poped.db)
FIM <- output$ret
det(FIM)
