# Section 1.2 R commands

# Section 1.2.2

library(LearnBayes)
data(studentdata)
studentdata[1,]
attach(studentdata)

# Section 1.2.3

table(Drink)
barplot(table(Drink),xlab="Drink",ylab="Count")

S=readline(prompt="Type  <Return>   to continue : ")

windows()
hours.of.sleep = WakeUp - ToSleep
summary(hours.of.sleep)
hist(hours.of.sleep,main="")

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.2.4

windows()

boxplot(hours.of.sleep~Gender,
  ylab="Hours of Sleep")

female.Haircut=Haircut[Gender=="female"]
male.Haircut=Haircut[Gender=="male"]
summary(female.Haircut)
summary(male.Haircut)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.2.5

windows()

plot(jitter(ToSleep),jitter(hours.of.sleep))

fit=lm(hours.of.sleep~ToSleep)
fit
abline(fit)
