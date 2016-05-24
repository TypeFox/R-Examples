
#for the FO approximation
ind=1

# no occasion defined in this example, so result is zero
output <- mf5(model_switch=t(poped.db$design$model_switch[ind,,drop=FALSE]),
   xt=t(poped.db$design$xt[ind,,drop=FALSE]),
   x=zeros(0,1),
   a=t(poped.db$design$a[ind,,drop=FALSE]),
   bpop=poped.db$parameters$bpop[,2,drop=FALSE],
   d=poped.db$parameters$param.pt.val$d,
   sigma=poped.db$parameters$sigma,
   docc=poped.db$parameters$param.pt.val$docc,
   poped.db)

# in this simple case the full FIM is just the sum of the individual FIMs
# and all the individual FIMs are the same
det(output$ret*32) == det(evaluate.fim(poped.db,fim.calc.type=4))  
