f1.glm2 <- 
    glm( factor(t2d) ~ Gdose + sex, fusion1m, family=binomial )
###hop:3-5
f1.glm2
###hop:3-9
summary(f1.glm2)
