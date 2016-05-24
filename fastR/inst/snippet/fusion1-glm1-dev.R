f1.glm0 <- glm( factor(t2d) ~ 1, fusion1m, family=binomial)
deviance(f1.glm0)
deviance(f1.glm1)
df1 <- df.residual(f1.glm0) - df.residual(f1.glm1); df1
1 - pchisq(deviance(f1.glm0) - deviance(f1.glm1), df=df1)
