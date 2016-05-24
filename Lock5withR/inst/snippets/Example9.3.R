head(RestaurantTips)
summary(lm(Tip ~ Bill, data = RestaurantTips)) 
confint(lm(Tip ~ Bill, data = RestaurantTips) , "Bill", level = 0.90)

