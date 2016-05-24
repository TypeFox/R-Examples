#####################################################
#  Section 7.2  Introduction to Hierarchical Modeling
#####################################################

library(LearnBayes)
library(lattice)

data(sluggerdata)

# fit logistic model for home run data for a particular player

logistic.fit=function(player)
{
d=subset(sluggerdata,Player==player)
x=d$Age; x2=d$Age^2
response=cbind(d$HR, d$AB-d$HR)
list(Age=x, p=glm(response~x+x2,family=binomial)$fitted)
}

names=unique(sluggerdata$Player); newdata=NULL
for (j in 1:9)
{
  fit=logistic.fit(as.character(names[j]))
  newdata=rbind(newdata,data.frame(as.character(names[j]),fit$Age,fit$p))
}
names(newdata)=c("Player","Age","Fitted")
xyplot(Fitted~Age|Player, data=newdata, type="l",lwd=3,col="black")

