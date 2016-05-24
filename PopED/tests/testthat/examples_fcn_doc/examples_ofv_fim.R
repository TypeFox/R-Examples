
## evaluate initial design 
FIM <- evaluate.fim(poped.db) # new name for function needed
FIM
get_rse(FIM,poped.db)

det(FIM)
ofv_fim(FIM,poped.db,ofv_calc_type=1) # det(FIM)
ofv_fim(FIM,poped.db,ofv_calc_type=2) # 1/trace_matrix(inv(FIM))
ofv_fim(FIM,poped.db,ofv_calc_type=4) # log(det(FIM)) 
ofv_fim(FIM,poped.db,ofv_calc_type=6) # Ds with fixed effects as "important"
ofv_fim(FIM,poped.db,ofv_calc_type=6,
        ds_index=c(1,1,1,0,0,0,1,1)) # Ds with random effects as "important"
ofv_fim(FIM,poped.db,ofv_calc_type=7) # 1/sum(get_rse(FIM,poped.db,use_percent=FALSE))

