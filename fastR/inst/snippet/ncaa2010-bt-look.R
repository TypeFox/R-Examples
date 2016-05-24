# look at top teams 
coef(ncaa.model)[rev(order(coef(ncaa.model)))[1:6]]
