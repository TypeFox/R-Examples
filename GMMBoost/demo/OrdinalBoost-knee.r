library(GMMBoost)
data(knee)

# sequential models

glm1 <- OrdinalBoost(pain ~ time + th + age + sex, rnd = list(id=~1),
        data = knee, model = "sequential", control = list(steps=10))
        
summary(glm1)

glm2 <- OrdinalBoost(pain ~ time + I(time^2) + th + age + sex, rnd = list(id=~1 + time),
        data = knee, model = "sequential", control = list(steps=10, katvar="age", lin="time",
        method="REML", sel.method="bic"))

summary(glm2)

# cumulative models
glm3 <- OrdinalBoost(pain ~ time + th + age + sex, rnd = list(id=~1), data = knee,
        model = "cumulative", control = list(steps=10))
        
summary(glm3)


# predict new data

knee.new<-knee[1:15,]
knee.new$id<-as.character(knee.new$id)
knee.new$id[1:3]<-c("1909","1909","1909")
knee.new$id<-as.factor(knee.new$id)
set.seed(123)
knee.new$age<-knee.new$age+rep(floor(rnorm(5,0,2)),each=3)

predict(glm1,knee.new)



