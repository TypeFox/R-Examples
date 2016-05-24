###################################################
# Section 4.4 A Bioassay Experiment
###################################################

library(LearnBayes)

 x = c(-0.86, -0.3, -0.05, 0.73)
 n = c(5, 5, 5, 5)
 y = c(0, 1, 3, 5)
 data = cbind(x, n, y)

 glmdata = cbind(y, n - y)
 results = glm(glmdata ~ x, family = binomial)
 summary(results)

# when x = -.7, median and 90th percentile of p are (.2,.4)
# when x = +.6, median and 90th percentile of p are (.8, .95)

a1.b1=beta.select(list(p=.5,x=.2),list(p=.9,x=.5))
a2.b2=beta.select(list(p=.5,x=.8),list(p=.9,x=.98))

prior=rbind(c(-0.7, 4.68, 1.12),
            c(0.6, 2.10, 0.74))
data.new=rbind(data, prior)

# plot prior #######################################

plot(c(-1,1),c(0,1),type="n",xlab="Dose",ylab="Prob(death)")
lines(-0.7*c(1,1),qbeta(c(.25,.75),a1.b1[1],a1.b1[2]),lwd=4)
lines(0.6*c(1,1),qbeta(c(.25,.75),a2.b2[1],a2.b2[2]),lwd=4)
points(c(-0.7,0.6),qbeta(.5,c(a1.b1[1],a2.b2[1]),c(a1.b1[2],a2.b2[2])),
    pch=19,cex=2)
text(-0.3,.2,"Beta(1.12, 3.56)")
text(.2,.8,"Beta(2.10, 0.74)")
response=rbind(a1.b1,a2.b2)
x=c(-0.7,0.6)
fit = glm(response ~ x, family = binomial)
curve(exp(fit$coef[1]+fit$coef[2]*x)/
     (1+exp(fit$coef[1]+fit$coef[2]*x)),add=T)

#######################################################
S=readline(prompt="Type  <Return>   to continue : ")

windows()
mycontour(logisticpost,c(-3,3,-1,9),data.new,
  xlab="beta0", ylab="beta1")

s=simcontour(logisticpost,c(-2,3,-1,11),data.new,1000)
points(s)

S=readline(prompt="Type  <Return>   to continue : ")

windows()
 plot(density(s$y),xlab="beta1")

S=readline(prompt="Type  <Return>   to continue : ")

 theta=-s$x/s$y
windows()
 hist(theta,xlab="LD-50",breaks=20)

 quantile(theta,c(.025,.975))
