###hop:3-9
summary(lm(rFlexibility~rStrength, punting))
punt.plot <- xyplot(rStrength~rFlexibility, punting)
# if all we want is the correlation coefficient, we can get it directly
r <- cor(punting$rStrength, punting$rFlexibility); r
r^2
