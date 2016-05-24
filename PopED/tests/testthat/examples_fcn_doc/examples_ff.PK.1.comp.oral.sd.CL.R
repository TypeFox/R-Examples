
##  create plot of model without variability 
plot_model_prediction(poped.db)

## evaluate initial design
FIM <- evaluate.fim(poped.db) 
FIM
det(FIM)
get_rse(FIM,poped.db)

