# build a variable that makes the model easier to describe
airp$x <- with( airp, ((loc==2) + 0.5*(loc==3)) )
model3 <- lm(pollution~ 1 + x, airp)
anova(model3,model)
