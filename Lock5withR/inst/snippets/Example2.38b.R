NutritionStudy60 = subset(NutritionStudy, Age>59)
xyplot(Alcohol ~ Calories, data = subset(NutritionStudy60, Alcohol<25))
cor(Alcohol ~ Calories, data = subset(NutritionStudy60, Alcohol<25))

