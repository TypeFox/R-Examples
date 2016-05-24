# select the version of this data set in the faraway package
data(orings,package="faraway")        
orings$failure <- orings$damage != 0   # convert to binary response
orings.model <- 
    glm(failure~temp,data=orings,family=binomial(link=logit))
###hop:3-9
summary(orings.model)
