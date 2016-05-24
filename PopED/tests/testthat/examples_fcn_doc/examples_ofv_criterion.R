
## evaluate initial design 
FIM <- evaluate.fim(poped.db) # new name for function needed
FIM
get_rse(FIM,poped.db)

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=1),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=1) # det(FIM)

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=2),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=2) 

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=4),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=4)

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=6),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=6)

ofv_criterion(ofv_fim(FIM,poped.db,ofv_calc_type=7),
              length(get_unfixed_params(poped.db)[["all"]]),
              poped.db,
              ofv_calc_type=7) 

