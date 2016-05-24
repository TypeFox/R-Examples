
# very few samples
poped.db$settings$ED_samp_size=10
ed_mftot(model_switch=poped.db$design$model_switch,
         groupsize=poped.db$design$groupsize,
         ni=poped.db$design$ni,
         xtoptn=poped.db$design$xt,
         xoptn=poped.db$design$x,
         aoptn=poped.db$design$a,
         bpopdescr=poped.db$parameters$bpop,
         ddescr=poped.db$parameters$d,
         covd=poped.db$parameters$covd,
         sigma=poped.db$parameters$sigma,
         docc=poped.db$parameters$docc, 
         poped.db)["ED_ofv"]

  
