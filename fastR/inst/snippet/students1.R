summary(students)
###hop:3-9
model <- lm(ACT~SAT,students); summary(model)
confint(model)
confint(lm( act ~ I(satm+satv), gpa))
