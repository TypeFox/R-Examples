
## @knitr 
beets <- c(41, 40, 41, 42, 44, 35, 41, 36, 47, 45)
no_beets <- c(51, 51, 50, 42, 40, 31, 43, 45)
c(xbar1=mean(beets), xbar2=mean(no_beets), 
  sd1=sd(beets), sd2=sd(no_beets))


## @knitr echo=FALSE
require(aplpack)


## @knitr stem.leaf.back
stem.leaf.backback(beets, no_beets, rule.line="Sturges") 


## @knitr dotchart-beets, echo=FALSE,  out.width=doublewide
library(lattice)
dotplot(ind ~ values, stack(list(Beets=beets, "No beets"=no_beets)), cex=3)
boxplot(beets, no_beets, names=c("beets", "no beets"), horizontal=TRUE)


## @knitr eval=FALSE
## boxplot(no_beets, beets, names=c("no beets", "beets"),horizontal=TRUE)


## @knitr 
speed <- michelson$Speed; expt <- michelson$Expt
fourth <- speed[expt == 4]
fifth <-  speed[expt == 5]


## @knitr 
d_4 <- density(fourth)
d_5 <- density(fifth)
xrange <- range(c(d_4$x, d_5$x))
yrange <- range(c(d_4$y, d_5$y))


## @knitr density-4th-5th, fig.keep="none", eval=FALSE
## plot(d_4, xlim=xrange, ylim=yrange, xlab="densities", main="")
## lines(d_5, lty=2)


## @knitr eval=FALSE
## legend(650, 0.008, legend=c("Fourth", "Fifth"), lty=c(1,2))


## @knitr michelson-qqplot, eval=FALSE
## qqplot(fourth, fifth)


## @knitr fig.keep="none"
ps <- seq(0.05, 0.95, by=0.05)
x <- quantile(fourth, ps)
y <- quantile(fifth, ps)
plot(x, y)                              # makes a scatter plot


## @knitr echo=FALSE,  out.width=doublewide
plot(d_4, xlim=xrange, ylim=yrange, xlab="densities", main="")
lines(d_5, lty=2)
qqplot(fourth, fifth)


## @knitr 
levels(ToothGrowth$supp)


## @knitr 
len_oj <- ToothGrowth$len[ToothGrowth$supp == "OJ"]
len_vc <- ToothGrowth$len[ToothGrowth$supp == "VC"]


## @knitr fig.keep="none"
boxplot(len_oj, len_vc, names=c("OJ", "VC"))


## @knitr 
first <- speed[expt == 1]
second <-  speed[expt == 1]
d_1 <- density(first)
d_2 <- density(second)
xrange <- range(c(d_1$x, d_2$x))
yrange <- range(c(d_1$y, d_2$y))
plot(d_1, xlim=xrange, ylim=yrange, xlab="densities", main="")
lines(d_2, lty=2)


## @knitr fig.keep="none"
qqplot(rivers, islands)                 # not straigh


## @knitr 
beets <- c(41, 40, 41, 42, 44, 35, 41, 36, 47, 45)
no_beets <- c(51, 51, 50, 42, 40, 31, 43, 45)
b <- list(beets = beets, "no beets" = no_beets)


## @knitr 
b                                       # print the list b


## @knitr 
b$beets


## @knitr 
b[1]                        # a list with just the first component
b[[1]]                                  # just the first component


## @knitr fig.keep="none", eval=FALSE
## boxplot(b)


## @knitr michelson-speeds, eval=FALSE
## speed <- michelson$Speed; expt <- michelson$Expt
## speeds <- list(speed[expt == 1], speed[expt == 2], speed[expt == 3],
##                speed[expt == 4], speed[expt == 5])
## names(speeds) <- paste("Expt", 1:5)
## boxplot(speeds)


## @knitr echo=FALSE, fig.height=4, fig.width=7, out.width=singlewide
speed <- michelson$Speed; expt <- michelson$Expt
speeds <- list(speed[expt == 1], speed[expt == 2], speed[expt == 3],
               speed[expt == 4], speed[expt == 5])
names(speeds) <- paste("Expt", 1:5)
boxplot(speeds)


## @knitr 
student <- c("Alice", "Bob")
grade <- c("A", "B")
attendance <- c("awesome", "bad")
data.frame(student, grade, attendance)


## @knitr fig.keep="none"
boxplot(Speed ~ Expt, data=michelson)


## @knitr fig.keep="none"
boxplot(Speed ~ Expt, data=michelson, subset=Run %in% 11:20)


## @knitr fig.keep="none"
plot(Speed ~ Expt, data=michelson)      # boxplots: Expt is a factor


## @knitr plot-summary-formula, eval=FALSE
## out <- summary(Speed ~ Expt, data=michelson)
## plot(out)


## @knitr echo=FALSE, out.width=singlewide
out <- summary(Speed ~ Expt, data=michelson) 
plot(out)


## @knitr 
b <- list("beets" = beets, "no beets" = no_beets)
stacked <- stack(b)
str(stacked)                            # two-variable data frame
headtail(stacked)                       # visual of data frame


## @knitr fig.keep="none"
plot(values ~ ind, data=stacked)


## @knitr 
speed <- michelson$Speed; expt <- michelson$Expt
speeds <- list(speed[expt == 1], speed[expt == 2], speed[expt == 3],
               speed[expt == 4], speed[expt == 5])
names(speeds) <- paste("Expt", 1:5, sep=" ")


## @knitr 
speeds <- split(michelson$Speed, michelson$Expt)
names(speeds) <- paste("Expt", 1:5, sep=" ")


## @knitr fig.keep="none"
HW_time <- list(Marsha= c(25, 0, 45, 90, 0), 
                Bill=c(30, 30, 30, 30), 
                Holly=c(15, 0, 90, 0))
boxplot(HW_time)


## @knitr 
vals <- data.frame(temp=c(98.2, 98.6, 98.2),
                   pulse=c(96, 56, 76),
                   systolic=c(134, 120, 150),
                   diastolic=c(90, 80, 95))
rownames(vals) <- c("Jackie", "Florence", "Mildred")


## @knitr eval=FALSE
## require(Hmisc)                          # loaded by UsingR
## library(MASS)
## summary(Speed ~ Expt, michelson)


## @knitr 
split(mtcars$mpg, mtcars$cyl)


## @knitr 
split(ToothGrowth$len, list(ToothGrowth$supp, ToothGrowth$dose))


## @knitr fig.keep="none"
time <- reaction.time$time; control <- reaction.time$control
boxplot(time[control == "C"], time[control == "T"])


## @knitr fig.keep="none"
boxplot(twins$Foster, twins$Biological, names=c("Foster","Biological"))


## @knitr 
plot(density(stud.recs$sat.m))
lines(density(stud.recs$sat.v),lty=2)
qqplot(stud.recs$sat.m, stud.recs$sat.v)


## @knitr 
temperature <- normtemp$temperature; gender <- normtemp$gender
male <- temperature[gender == 1]
female <- temperature[gender == 2]
boxplot(list(male=male, female=female))


## @knitr plot-wrist-height, fig.keep="none"
plot(height ~ wrist, data=fat, subset=height > 50)


## @knitr plotoptions-table, results="asis", echo=FALSE
tbl <- "
\\code{main}: Title to put on the graphic.
\\code{xlab}: Label for the $x$-axis. Similarly for \\code{ylab}.
\\code{xlim}: Specify the $x$-limits, as in \\code{xlim=c(0,10)}, for the interval $[0,10]$. imilarly for \\code{ylim}.
\\code{type}: Type of plot to make. Use \\code{\"p\"} for
              points (the default), \\code{\"l\"} (ell not one) for lines, and
              \\code{\"h\"} for vertical lines.
\\code{bty}:  Type of box to draw. Use \\code{\"l\"} for
              ``L''-shaped, default is \\code{\"o\"}, which is ``O''-shaped.
             Details in \\code{?par}. 
\\code{pch}: Plot character used. This can be a number or a single character. Numbers $0$-$18$ 
             are similar to those of \\proglang{S}, though there are more available beyond that range.
\\code{lty}: When lines are plotted, specifies the type of line
             to be drawn. Details in \\code{?par}. 
\\code{lwd}: The thickness of lines. Numbers bigger than 1 increase the default.
\\code{cex}: Magnification factor. Default is 1.
\\code{col}: Specifies the color to use for the points or lines.  
"
out <- read.dcf(textConnection(tbl))
d <- data.frame(Command=colnames(out), Description=out[1,], stringsAsFactors=FALSE)
cat(booktabs(d, 
         caption="Useful arguments for various plotting commands.", 
         subtop=NULL, 
         rnames=NULL, 
         label="tab:plot-arguments",
         colTypes=c("l", "p{3.5in}")
         ))


## @knitr 
names(fat)


## @knitr 
mean(fat$neck) / mean(fat$wrist)        # ratio of means
mean(fat$neck/fat$wrist)                # mean of ratios


## @knitr plot-neck-wrist-xy, fig.keep="none"
plot(fat$wrist, fat$neck)


## @knitr plot-neck-wrist, fig.keep="none"
plot(neck ~ wrist, data=fat)


## @knitr plot-neck-wrist-twenties, fig.keep="none"
plot(neck ~ wrist, data=fat, subset = 20 <= age &  age < 30)


## @knitr echo=FALSE, out.width=doublewide
plot(neck ~ wrist, data=fat)
plot(height ~ wrist, data=fat, subset=height > 50)


## @knitr crosshairs, fig.keep="none"
x <- fat$wrist[1:20]; y <- fat$neck[1:20] # some data
plot(x, y, main="Neck by wrist")
abline(v = mean(x), lty=2)                # dashed vertical line
abline(h = mean(y), lty=2)                # dashed horizontal line
points(mean(x), mean(y), pch=16, cex=4, col=rgb(0,0,0, .25))


## @knitr plot-layers-table, results="asis", echo=FALSE
tbl <- "
\\code{points}: Add points to a graphic.
\\code{lines}: Add line segments to a graphic.
\\code{abline}: Add a line to a graphic.
\\code{text}: Locate text in body of a graphic.
\\code{mtext}: Locate text in margins of a graphic.
\\code{title}: Add main title, $x$ label or $y$ label to a graphic
\\code{legend}: Locate legend in body of a graphic.
"
out <- read.dcf(textConnection(tbl))
d <- data.frame(Function=colnames(out), Description=out[1,], stringsAsFactors=FALSE)
cat(booktabs(d, 
         caption="Plotting functions to add to an existing graphic.", 
         subtop=NULL, 
         rnames=NULL, 
         label="tab:plot-layers",
         colTypes=c("l", "p{3.5in}")
         ))


## @knitr echo=FALSE, out.width=doublewide
x <- fat$wrist[1:20]; y <- fat$neck[1:20] # some data
plot(x, y, main="Neck by wrist")
abline(v = mean(x), lty=2)                # dashed vertical line
abline(h = mean(y), lty=2)                # dashed horizontal line
points(mean(x), mean(y), pch=16, cex=4, col=rgb(0,0,0, .25))
x <- fat$age[1:20]; y <- fat$ankle[1:20]
plot(x, y, main="Ankle by age")
abline(v = mean(x), lty=2)                # dashed vertical line
abline(h = mean(y), lty=2)                # dashed horizontal line
points(mean(x), mean(y), pch=16, cex=4, col=rgb(0,0,0, .25))


## @knitr 
cor(fat$wrist, fat$neck)                # strongly correlated
cor(fat$wrist, fat$height)              # mildly correlated
cor(fat$age, fat$ankle)                 # basically uncorrelated


## @knitr 
cor(Animals$body, Animals$brain)


## @knitr 
body <- Animals$body; brain <- Animals$brain 
cross_prods <- (body - mean(body)) * (brain - mean(brain))
Animals[cross_prods < 0, ]


## @knitr 
cor(rank(body), rank(brain))


## @knitr 
cor(body, brain, method="spearman")


## @knitr 
cor(ToothGrowth$dose, ToothGrowth$len)


## @knitr 
l <- split(ToothGrowth$len, ToothGrowth$dose)
group_means <- c(mean(l[[1]]), mean(l[[2]]), mean(l[[3]]))
cor(c(0.5, 1, 2), group_means)


## @knitr 
cor(SAT$salary, SAT$total)


## @knitr SAT-mark-perc, eval=FALSE
## plot(  total ~ salary, SAT)
## points(total ~ salary, SAT, subset = perc < 10, pch=15)  ## square
## points(total ~ salary, SAT, subset = perc > 40, pch=16)  ## solid


## @knitr echo=FALSE, out.width=singlewide
plot(total ~ salary, SAT)
plot(  total ~ salary, SAT)
points(total ~ salary, SAT, subset = perc < 10, pch=15)  ## square
points(total ~ salary, SAT, subset = perc > 40, pch=16)  ## solid


## @knitr 
total <- SAT$total; salary <- SAT$salary; perc <- SAT$perc
less_10 <- perc < 10
more_40 <- perc > 40
between <- !less_10 & !more_40
c(less = cor(total[less_10], salary[less_10]),
  between = cor(total[between], salary[between]),
  more = cor(total[more_40], salary[more_40]))


## @knitr 
buyback <- c(ACT=1500, NSW=2500, WA=2700, Qld=3700, Vic=4250, 
             SA=4200, NT=5000, Tas=7500)
change <- -c(.2, 1,.9, 2.5, 1.0, 1.6, 2.5, 3.4)
m <- data.frame(buyback, change)
m


## @knitr 
cor(buyback, change)


## @knitr ToothGrowth-mean-line, fig.keep="none"
plot(len ~ dose, data=ToothGrowth, pch=16, col=rgb(0,0,0,.25))
points(c(0.5, 1, 2), group_means, cex=1.5, pch=18)
lines(c(0.5, 1, 2), group_means)


## @knitr echo=FALSE, out.width=doublewide
plot(len ~ dose, data=ToothGrowth, pch=16, col=rgb(0,0,0,.25))
points(c(0.5, 1, 2), group_means, cex=1.5, pch=18)
lines(c(0.5, 1, 2), group_means)


## @knitr echo=FALSE, out.width=singlewide
x <- c(1,  2,  3, 3, 4,  5)
y <- c(3, -2, -1, 2, 2, 5)

plot(y ~ x, bty="l", pch=16)
res <- lm(y ~ x)
abline(res)

bs <- coef(res)

yhat <- predict(res, data.frame(x=3))
## annotate
points(3, yhat, pch=18, cex=3, col=rgb(0,0,0, .5))

text(3 - .1, -1   , expression(paste("(", x[3], ", ", y[3],")", sep="")),     pos=2)
text(3 - .1, yhat , expression(paste("(", x[3], ", ", hat(y)[3], ")", sep="")),pos=2)

lines(c(3,3), c(-1 + .1, yhat - .1), lty=2)

text(3 + .1, mean(c(-1, yhat)), expression(residual == y[3] - hat(y)[3]), pos=4)




## @knitr 
res <- lm(maxrate ~ age, data=heartrate)
res


## @knitr fig.keep="none"
plot(maxrate ~ age, data=heartrate)
abline(res)                             # add line to graphic


## @knitr 
sum(resid(res))


## @knitr 
actual <- heartrate$maxrate
diffs <- resid(res) - (actual - fitted(res))
max(diffs)


## @knitr 
age <- c(30, 40)
coef(res)[1] + coef(res)[2] * age       # beta_0 + beta_1 x


## @knitr 
predict(res, data.frame(age=c(30, 40)))


## @knitr extractor-functions-table, eval=FALSE, echo=FALSE
## tbl <- "
## \\code{coef}: Returns estimated coefficients.
## \\code{residuals}: Returns residual (observed $-$ expected).
## \\code{fitted}: Returns fitted response values ($\\HAT{y})$.
## \\code{predict}: Predict fitted value for new inputs.
## \\code{plot}: Produce plot(s) to assess model assumptions.
## \\code{summary}: Summarize model fit for inference.
## \\code{anova}: Comparison of models.
## \\code{model.matrix}: Matrix used to fit model mathematically.
## \\code{vcov}: Estimated covariance between parameter estimates.
## "
## out <- read.dcf(textConnection(tbl))
## d <- data.frame(Function=colnames(out), Description=out[1,], stringsAsFactors=FALSE)
## cat(booktabs(d,
##          caption="Non-exhaustive listing of extractor functions for working with R's modeling objects.",
##          subtop=NULL,
##          rnames=NULL,
##          label="tab:extractor-functions",
##          colTypes=c("l", "p{3.5in}")
##          ))


## @knitr brain-body-log-log, eval=FALSE
## plot(log(brain) ~ log(body), data=Animals)
## ## label the outliers
## idx <-  log(Animals$body) > 6 & log(Animals$brain) < 6 # select
## text(log(Animals$body[idx]), log(Animals$brain[idx]), # label
##      labels=rownames(Animals)[idx], pos=2)


## @knitr echo=FALSE
## not evaled in plot above, do so here
idx <-  log(Animals$body) > 6 & log(Animals$brain) < 6 # select


## @knitr cache=FALSE
lm(log(brain) ~ log(body), data=Animals, subset=!idx)


## @knitr echo=FALSE, out.width=singlewide
plot(log(brain) ~ log(body), data=Animals)
## label the outliers
idx <-  log(Animals$body) > 6 & log(Animals$brain) < 6 # select
text(log(Animals$body[idx]), log(Animals$brain[idx]), # label
     labels=rownames(Animals)[idx], pos=2)


## @knitr 
idx <- 60 <= kid.weights$age & kid.weights$age < 72


## @knitr fig.out="none"
plot(weight/2.2 ~ (height*2.54/100)^2, data=kid.weights, subset=idx)


## @knitr 
fm <- weight/2.2 ~ I(height*2.54/100^2)


## @knitr fig.keep="none"
plot(fm, data=kid.weights, subset=idx)
res <- lm(fm, data=kid.weights, subset=idx)
abline(res)
res


## @knitr echo=FALSE, out.width=singlewide
parent <- father.son$fheight
child <- father.son$sheight
plot(child ~ parent, col=rgb(0,0,0, alpha=0.15), pch=16, main="Pearson's data")
abline(lm(child ~ parent), lty=1, lwd=2)
## use SD(child)/SD(parent) for slope
curve(mean(child) + sd(child)/sd(parent) * (x - mean(parent)), add=TRUE, lty=2, lwd=2)
legend(71.5, 64, c("lm", "sd"), lty=1:2, lwd=2)



## @knitr quantreg-engel, eval=FALSE
## require(quantreg)                       # also loaded with UsingR
## data(engel)                             # not loaded by default
## f <- foodexp ~ income
## res.lm <- lm(f,          data=engel)
## res.rq <- rq(f, tau=0.5, data=engel)
## plot(f, data=engel)
## abline(res.lm, lty=1)
## abline(res.rq, lty=2)
## legend(4000, 1500, c("lm", "rq"), lty=1:2)


## @knitr engel_data, echo=FALSE, out.width=doublewide
require(quantreg)                       # also loaded with UsingR
data(engel)                             # not loaded by default
f <- foodexp ~ income
res.lm <- lm(f,          data=engel)
res.rq <- rq(f, tau=0.5, data=engel)
plot(f, data=engel)
abline(res.lm, lty=1)
abline(res.rq, lty=2)
legend(4000, 1500, c("lm", "rq"), lty=1:2)
e <- transform(subset(engel, income <= 3000), income_level = cut(income, seq(0, 3000, by=250)))
p <- ggplot(e, aes(y=foodexp, x=income, group=income_level)) + geom_boxplot() + geom_quantile(aes(group=NULL), quantiles=0.5)
suppressWarnings(suppressMessages(print(p)))


## @knitr echo=FALSE
## Function to minimize distance to a line
## y = b0 + b1 x --> 0 = b1 x - 1 y + b0
## 

d2 <- function(b0,b1,x,y) {
  a = b1; b = -1; c=b0
  abs(a * x  + b*y +c)/sqrt(a^2 + b^2)
}
minimize_it <- function(a,b) {
  g = function(b0, b1) {
    sum(mapply(function(x,y) d2(b0, b1, x, y), a, b))
  }
  fr = function(x) g(x[1], x[2])

  out = optim(c(0,0), fr)
  
  out$par
}
  


## @knitr echo=FALSE, out.width=singlewide
parent <- father.son$fheight
child <- father.son$sheight
plot(child ~ parent, col=rgb(0,0,0, alpha=0.15), pch=16)
abline(lm(child ~ parent), lty=1, lwd=2)
## use SD(child)/SD(parent) for slope
#curve(mean(child) + sd(child)/sd(parent) * (x - mean(parent)), add=TRUE, lty=2, lwd=2)
## princomp
#out <- princomp(~ child + parent)
#m <- Reduce("/", out$loadings[,1])
out <- prcomp(~ child + parent)
m <- Reduce("/",out$rotation[,1])
curve(mean(child) + m * (x - mean(parent)), add=TRUE, lty=2, lwd=2)
legend(71.5, 64, c("lm", "PCA"), lty=1:2, lwd=2)


## @knitr animals-rlm-lm, eval=FALSE
## f <- log(brain) ~ log(body)
## res.lm <- lm(f, data=Animals)
## res.rlm <- rlm(f, data=Animals)
## plot(f, data=Animals)
## abline(res.lm, lty=1)
## abline(res.rlm, lty=2)
## legend(-2, 8, legend=c("lm", "rlm"), lty=1:2)


## @knitr echo=FALSE, out.width=singlewide
f <- log(brain) ~ log(body)
res.lm <- lm(f, data=Animals)
res.rlm <- rlm(f, data=Animals)
plot(f, data=Animals)
abline(res.lm, lty=1)
abline(res.rlm, lty=2)
legend(-2, 8, legend=c("lm", "rlm"), lty=1:2)


## @knitr 
which(resid(res.lm) < -2)


## @knitr loess-fit, eval=FALSE
## f <- child ~ parent
## plot(f, data=Galton, col=rgb(0,0,0, alpha=0.25))
## ##
## res.lm <- lm(f, data=Galton)
## abline(res.lm, lty=1, lwd=2)
## ##
## res.loess <- loess(f, data=Galton, degree=1)  # line, not quadratic
## rng <- seq(64, 73, length=20)
## newdata <- data.frame(parent=rng)
## predicted <- predict(res.loess, newdata)      # predicted  by loess
## lines(predicted ~ rng, lty=2, lwd=2)


## @knitr echo=FALSE, out.width=singlewide
f <- child ~ parent
plot(f, data=Galton, col=rgb(0,0,0, alpha=0.25))
##
res.lm <- lm(f, data=Galton)
abline(res.lm, lty=1, lwd=2)
##
res.loess <- loess(f, data=Galton, degree=1)  # line, not quadratic
rng <- seq(64, 73, length=20)
newdata <- data.frame(parent=rng)
predicted <- predict(res.loess, newdata)      # predicted  by loess
lines(predicted ~ rng, lty=2, lwd=2)


## @knitr fig.keep="none"
hd_ratio <- homedata$y2000 / homedata$y1970
hist(hd_ratio, probability=TRUE)
lines(density(hd_ratio))


## @knitr 
plot(hr ~ temperature, normtemp)
cor(normtemp$temperature, normtemp$hr)


## @knitr 
plot(Biological ~ Foster, data=twins)


## @knitr 
cor(twins$Foster, twins$Biological)
cor(twins$Foster, twins$Biological, method="spearman")


## @knitr 
x77 <- data.frame(state.x77);


## @knitr 
plot(Population ~ Frost, data=x77)
plot(Murder ~ Population, data=x77)
plot(Population ~ Area, data=x77)
plot(Income ~ HS.Grad, data=x77)      


## @knitr 
names(nym.2002)                         # variable names
cor(nym.2002$age,  nym.2002$time)


## @knitr fig.keep="none"
plot(y ~ x, data=state.center)


## @knitr 
cor(batting$HR, batting$SO)


## @knitr 
res <- lm(abdomen ~ wrist, data = fat)
predict(res, data.frame(wrist=17))


## @knitr fig.keep="none"
cor(wtloss$Weight, wtloss$Days)             # close to -1
plot(Weight ~ Days, data=wtloss)          
res <- lm(Weight ~ Days, data=wtloss)
abline(res)                   
plot(residuals(res) ~ Days, data=wtloss)    # A trend


## @knitr 
x77 <- data.frame(state.x77)    


## @knitr 
plot(Illiteracy ~ HS.Grad, data = x77)      
plot(Life.Exp ~ Murder, data = x77)
plot(Income ~ Illiteracy, data = x77)


## @knitr 
res <- lm(RBI ~ HR, data=batting)
predict(res, data.frame(HR=33))
98 - predict(res, data.frame(HR=33)) 


## @knitr fig.keep="none"
lm(Female ~ Male, data = too.young)
plot(Female ~ Male, data = too.young)
abline(lm(Female ~ Male, data = too.young))
abline(7,1/2,lwd=2)


## @knitr fig.keep="none"
plot(price ~ carat, data = diamond, pch=5)    
res <- lm(price ~ carat, data = diamond)
abline(res)
predict(res, data.frame(carat=1/3))


## @knitr 
plot(log(time) ~ voltage, data=breakdown)
res <- lm(log(time) ~ voltage, data=breakdown)
abline(res)
res


## @knitr fig.keep="none"
data(motors, package="MASS")
plot(time ~ temp, pch=cens, data=motors)


## @knitr 
table(motors$temp)


## @knitr 
res <- lm(time ~ temp,  data=motors)
res


## @knitr 
predict(res, newdata=data.frame(temp=210))


## @knitr 
rbind(c(56,8), c(2,16))        # combine rows
cbind(c(56,2), c(8,16))        # bind as columns


## @knitr 
seatbelts <- matrix(c(56, 2, 8, 16), nrow=2)


## @knitr eval=FALSE
## x <- matrix(1)                 # need to initialize x
## x <- edit(x)                   # will edit matrix with spreadsheet


## @knitr 
rownames(seatbelts) <- c("buckled","unbuckled")
colnames(seatbelts) <- c("buckled","unbuckled")
seatbelts


## @knitr 
dimnames(seatbelts) <- list(parent=c("buckled","unbuckled"), 
                            child=c("buckled","unbuckled")) 


## @knitr 
x <- c(56,8); names(x) = c("buckled","unbuckled")
y <- c(2,16)
rbind(buckled=x, unbuckled=y)  # names rows, columns come from x


## @knitr 
headtail(grades)
table(grades$prev, grades$grade)        # also: table(grades)


## @knitr 
seatbelts
margin.table(seatbelts, margin=1)       # row sum is for parents
margin.table(seatbelts, margin=2)       # column sum for kids


## @knitr 
addmargins(seatbelts)


## @knitr 
tbl <- with(grades, table(prev, grade))
margin.table(tbl, 1) 
margin.table(tbl, 2)


## @knitr echo=FALSE
options(digits=1)


## @knitr 
prop.table(table(grades$prev, grades$grade), margin=1) * 100


## @knitr echo=FALSE
options(digits=4)


## @knitr 
headtail(Fingerprints)


## @knitr 
idx <- !is.na(Fingerprints$count)       # issue with NA in rep
whorls <- rep(Fingerprints$Whorls[idx], Fingerprints$count[idx])
loops <- rep(Fingerprints$Loops[idx], Fingerprints$count[idx])
table(whorls, loops)


## @knitr 
xtabs(count ~ Whorls + Loops, Fingerprints)


## @knitr 
xtabs( ~ Origin + Type, Cars93)          # no LHS, so tallies


## @knitr 
tbl <- xtabs( ~ Origin + Type +  Passengers , Cars93)
ftable(tbl, row.vars=2, col.vars=c(1,3))


## @knitr 
price.range <- cut(Cars93$Price, c(0,  20, 70))
tbl <- xtabs( ~ price.range + Origin + Cylinders, data=Cars93)
out <- ftable(tbl, row.vars=3)
100 * prop.table(out, margin=1)


## @knitr seatbelt-barplots, eval=FALSE
## barplot(seatbelts, xlab="Parent", main="Child seat-belt usage")
## barplot(seatbelts, xlab="Parent", main="Child seat-belt usage",
##         beside=TRUE)


## @knitr echo=FALSE, out.width=doublewide
barplot(seatbelts, xlab="Parent", main="Child seat-belt usage")
barplot(seatbelts, xlab="Parent", main="Child seat-belt usage",
        beside=TRUE)


## @knitr 
titanic <- as.data.frame(Titanic)
xtabs(Freq ~ Survived + Class, data=titanic, subset=Sex=="Female")


## @knitr mosaicplot-titanic, eval=FALSE
## tbl <- xtabs(Freq ~ Sex, titanic)
## mosaicplot(tbl)
## tbl <- xtabs(Freq ~ Sex + Survived, titanic)
## mosaicplot(tbl)


## @knitr echo=FALSE, out.width=doublewide
tbl <- xtabs(Freq ~ Sex, titanic)
mosaicplot(tbl)
tbl <- xtabs(Freq ~ Sex + Survived, titanic)
mosaicplot(tbl)


## @knitr mosaicplot-titanic-with-class, eval=FALSE
## tbl <- xtabs(Freq ~ Sex + Survived  + Class, titanic)
## mosaicplot(tbl)


## @knitr fig.keep="none"
mosaicplot( xtabs(Freq ~ Class + Sex + Survived, data=titanic))


## @knitr echo=FALSE, out.width=singlewide
tbl <- xtabs(Freq ~ Sex + Survived  + Class, titanic)
mosaicplot(tbl)


## @knitr titanic-survive-class, echo=FALSE, eval=FALSE
## tbl <- xtabs(Freq ~ Survived + Class, titanic)
## mosaicplot(tbl)


## @knitr echo=FALSE, out.width=singlewide
tbl <- xtabs(Freq ~ Survived + Class, titanic)
mosaicplot(tbl)


## @knitr 
Survived <- rep(titanic$Survived, titanic$Freq)
Survived <- ordered(Survived)
Class <-  rep(titanic$Class, titanic$Freq)
Class <- ordered(Class)
head(Class)


## @knitr 
cor(as.numeric(Survived), as.numeric(Class), method="kendall")


## @knitr kendall, echo=FALSE, out.width=singlewide
set.seed(2.71)
ox <- rnorm(5); oy <- rnorm(5)
x <- rank(ox); y <- rank(oy)
plot(x,y, pch=16, cex=2)
i <- 2
abline(v=x[i], lty=2)
abline(h=y[i], lty=2)
concordant <- sum((x[-i]-x[i])*(y[-i] - y[i]) > 0)
discordant <- sum((x[-i]-x[i])*(y[-i] - y[i]) < 0)
title(paste("Concordant=", concordant, ", discordant=", discordant, ", tau=", cor(ox,oy,method="kendall"))) 


## @knitr 
f <- Freq ~ Survived + Class 
tbl <- xtabs(f, data=titanic, subset=Sex=="Female")
summary(tbl)


## @knitr 
margin.table(tbl,1)
ptop <- 126/sum(margin.table(tbl,1))
margin.table(tbl, 2)
pleft <- 145/sum(margin.table(tbl,2)) 
fe <- ptop * pleft * sum(tbl)
fo <- tbl[1,1]
c(fo=fo, fe=fe, residual=fo - fe)       # print


## @knitr 
chisq.test(tbl)$expected


## @knitr 
fo <- tbl
fe <- chisq.test(tbl)$expected
(fo - fe)^2 / fe
sum((fo - fe)^2 / fe)


## @knitr 
2 * sum(fo * log(fo/fe))                # no fe can be 0


## @knitr 
xtabs(Freq ~ Sex + Survived, titanic)


## @knitr 
table(coins$value)            
table(coins$value) * c(0.01, 0.05, 0.1, 0.25) 
sum(table(coins$value) * c(0.01, 0.05, 0.1, 0.25))


## @knitr fig.keep="none"
barplot(table(coins$year))


## @knitr 
table(cut(coins$year, breaks = seq(1920,2010,by=10)))


## @knitr 
cor(coins$year, coins$value)


## @knitr 
require(MASS) 
table(UScereal$shelf, UScereal$mfr)


