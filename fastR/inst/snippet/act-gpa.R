grades <- actgpa
t.test(grades$ACT)
t.test(grades$GPA)
grades.model <- lm(GPA~ACT,grades)
summary(grades.model)
grades.plot1 <- xyplot(GPA~ACT,data=grades,panel=panel.lmbands)
predict(grades.model,new=data.frame(ACT=25),interval="confidence")
predict(grades.model,new=data.frame(ACT=30),interval="prediction")
