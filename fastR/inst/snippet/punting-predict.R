newdata <- data.frame(rStrength=175,rFlexibility=100)
predict(punting.lm, new=newdata, interval='confidence')
predict(punting.lm, new=newdata, interval='prediction')
