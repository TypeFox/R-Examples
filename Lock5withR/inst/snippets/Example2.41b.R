Tip.Fun <- makeFun(lm(Tip ~ Bill, data = RestaurantTips)) # make a function of the linear model
Tip.Fun(Bill = 59.33) # predicted tip when bill is $59.33
Tip.Fun(Bill = 9.52)
Tip.Fun(Bill = 23.70)

