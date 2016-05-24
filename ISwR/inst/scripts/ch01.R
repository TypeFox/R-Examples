plot(rnorm(500))
2 + 2
exp(-2)
rnorm(15)
x <- 2
x
x + x
weight <- c(60, 72, 57, 90, 95, 72)
weight
height <- c(1.75, 1.80, 1.65, 1.90, 1.74, 1.91)
bmi <- weight/height^2
bmi
sum(weight)
sum(weight)/length(weight)
xbar <- sum(weight)/length(weight)
weight - xbar
(weight - xbar)^2
sum((weight - xbar)^2)
sqrt(sum((weight - xbar)^2)/(length(weight) - 1))
mean(weight)
sd(weight)
t.test(bmi, mu=22.5)
plot(height,weight)
if (.make.epsf) dev.copy2eps(file="h-w.ps")
plot(height, weight, pch=2)
if (.make.epsf) dev.copy2eps(file="h-w-triangle.ps")
hh <- c(1.65, 1.70, 1.75, 1.80, 1.85, 1.90)
lines(hh, 22.5 * hh^2)
if (.make.epsf) dev.copy2eps(file="h-w-line.ps")
args(plot.default)
c("Huey","Dewey","Louie")
c('Huey','Dewey','Louie')
c(T,T,F,T)
bmi > 25
cat(c("Huey","Dewey","Louie"))
cat("Huey","Dewey","Louie", "\n")
cat("What is \"R\"?\n")
c(42,57,12,39,1,3,4)
x <- c(1, 2, 3)
y <- c(10, 20)
c(x, y, 5)
x <- c(red="Huey", blue="Dewey", green="Louie")
x
names(x)
c(FALSE, 3)
c(pi, "abc")
c(FALSE, "abc")
seq(4,9)
seq(4,10,2)
4:9
oops <- c(7,9,13)
rep(oops,3)
rep(oops,1:3)
rep(1:2,c(10,15))
x <- 1:12
dim(x) <- c(3,4)
x
matrix(1:12,nrow=3,byrow=T)
x <- matrix(1:12,nrow=3,byrow=T)
rownames(x) <- LETTERS[1:3]
x
t(x)
cbind(A=1:4,B=5:8,C=9:12)
rbind(A=1:4,B=5:8,C=9:12)
pain <- c(0,3,2,2,1)
fpain <- factor(pain,levels=0:3)
levels(fpain) <- c("none","mild","medium","severe")
fpain
as.numeric(fpain)
levels(fpain)
intake.pre <- c(5260,5470,5640,6180,6390,
6515,6805,7515,7515,8230,8770)
intake.post <- c(3910,4220,3885,5160,5645,
4680,5265,5975,6790,6900,7335)
mylist <- list(before=intake.pre,after=intake.post)
mylist
mylist$before
d <- data.frame(intake.pre,intake.post)
d
d$intake.pre
intake.pre[5]
intake.pre[c(3,5,7)]
v <- c(3,5,7)
intake.pre[v]
intake.pre[1:5]
intake.pre[-c(3,5,7)]
intake.post[intake.pre > 7000]
intake.post[intake.pre > 7000 & intake.pre <= 8000]
intake.pre > 7000 & intake.pre <= 8000
d <- data.frame(intake.pre,intake.post)
d[5,1]
d[5,]
d[d$intake.pre>7000,]
sel <- d$intake.pre>7000
sel
d[sel,]
d[1:2,]
head(d)
energy
exp.lean <- energy$expend[energy$stature=="lean"]
exp.obese <- energy$expend[energy$stature=="obese"]
l <- split(energy$expend, energy$stature)
l
lapply(thuesen, mean, na.rm=T)
sapply(thuesen, mean, na.rm=T)
replicate(10,mean(rexp(20)))
m <- matrix(rnorm(12),4)
m
apply(m, 2, min)
tapply(energy$expend, energy$stature, median)
intake$post
sort(intake$post)
order(intake$post)
o <- order(intake$post)
intake$post[o]
intake$pre[o]
intake.sorted <- intake[o,] 
save.image("ch1.RData") 
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
