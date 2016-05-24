coag <- coagulation
coag$x1 <- coag$diet=='B'
coag$x2 <- coag$diet=='C'
coag$x3 <- coag$diet=='D'
coag.model <- lm(coag~x1+x2+x3,coag)
coag.model1 <- lm(coag~1,coag)
anova(coag.model1,coag.model)
