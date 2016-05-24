# build two variables that make the model easier to describe
airp$x1 <- with( airp, ((loc==1) + 0.5*(loc==3)) )
airp$x2 <- with( airp, ((loc==2) + 0.5*(loc==3)) )
model3 <- lm(pollution~ -1 + x1 + x2, airp)
anova(model3,model)
