## print a shortend summary for lm model
short_summary <- function(x) {
  x <- summary(x)
  cmat <- coef(x)
  printCoefmat(cmat)
}



## @knitr echo=FALSE, out.width=singlewide
## make dotplot.R

x = c(2166, 1568, 2233, 1882, 2019)
y = c(2279, 2075, 2131, 2009, 1793)
z = c(2226 ,2154 ,2583 ,2010 ,2190)
values = c(x,y,z)
ind =  gl(3, 5, 15, labels=c("May", "Sep", "Dec"))
## draw chart
stripchart(values ~ ind,pch=16,cex=2,ylim=c(.7,3.2))
## label
arrows(mean(x),.9,mean(x),.91,lwd=2,length=.15)
text(mean(x),.7,expression(bar(May)),cex=1.5)

arrows(mean(y),1.9,mean(y),1.91,lwd=2,length=.15)
text(mean(y),1.7,expression(bar(Sep)),cex=1.5)
  
arrows(mean(z),2.9,mean(z),2.91,lwd=2,length=.15)
text(mean(z),2.7,expression(bar(Dec)),cex=1.5)

lines(c(80,110),c(1,1),lty=2)
lines(c(80,110),c(2,2),lty=2)
lines(c(80,110),c(3,3),lty=2)

grand.mean = mean(c(x,y,z))
lines(c(grand.mean,grand.mean),c(.7,3.2),lty=3)

arrows(grand.mean,1.4,grand.mean,1.41,lwd=2,length=.15)
text(grand.mean,1.2,"grand mean",cex=1.5)


## @knitr 
may <- c(2166, 1568, 2233, 1882, 2019)
sep <- c(2279, 2075, 2131, 2009, 1793)
dec <- c(2226, 2154, 2583, 2010, 2190)
#
xbar <- mean(c(may, sep, dec))
SST <- (15-1) * var(c(may, sep, dec))   # (n-1) * var(.) is SS
SSE <- (5-1) * var(may) + (5-1) * var(sep) + (5-1) * var(dec)
SSTr <- 5 * ((mean(may) - xbar)^2 + (mean(sep) - xbar)^2 + 
            (mean(dec) - xbar)^2)
#
c(SST=SST, SSTr=SSTr, SSE=SSE)
#
F.obs = (SSTr/(3-1)) / (SSE/(15-3))
F.obs
pf(F.obs, df1=3 - 1, df2=15 - 3, lower.tail=FALSE)


## @knitr 
d <- stack(list(may=may, sep=sep, dec=dec)) # need names for list
names(d)                                    # two variables
oneway.test(values ~ ind, data=d, var.equal=TRUE)


## @knitr 
res <- aov(values ~ ind, data = d)
res    


## @knitr 
sqrt(586720/12)


## @knitr 
summary(res)


## @knitr fig.keep="none"
UBP <- c(168.2, 161.4, 163.2, 166.7, 173.0, 173.3, 
         160.1, 161.2, 166.8)
grip.type <- gl(3, 3, 9, labels=c("classic", "integrated", "modern"))
boxplot(UBP ~ grip.type, ylab="Power (watts)", 
        main="Effect of cross country grip")


## @knitr 
res <- aov(UBP ~ grip.type)
summary(res)


## @knitr echo=FALSE, out.width=singlewide
boxplot(UBP ~ grip.type, ylab="Power (watts)", main="Effect of cross country grip")


## @knitr eval=FALSE, echo=FALSE
## ## R code to download and create data set movie_data_2011
## library(XML)
## library(lubridate)
## 
## tpl = "http://www.the-numbers.com/box-office-chart/weekend/2011/%s/%s"
## 
## 
## get_weekend <- function(date) {
##   m <- gsub(" ", "0", sprintf("%2d", month(date)))
##   d <- gsub(" ", "0", sprintf("%2d", day(date)))
##   u <- sprintf(tpl, m, d)
## 
##   x <- readHTMLTable(u)[[4]]
##   ## tidy up tables
##   ## drop first column
##   x <- x[-1]
##   ## drop per theatre
##   x <- x[-8]
##   ## rename
##   nms <- c("Previous", "Movie", "Distributor", "Genre", "Gross", "Change", "Theaters", "TotalGross", "Days")
##   names(x) <- nms
##   ## make Previous numeric. Use NA for new
##   make_numeric <- function(x) {
##     if (x == "new")
##       return(NA)
##     else if (x == "(-)")
##       return(NA)
##     else
##       as.numeric(gsub("\\(", "", gsub(")", "", as.character(x))))
##   }
##   x$Previous <- sapply(x$Previous, make_numeric)
## 
##   ## strip off $ and ,
##   make_numeric <- function(x) {
##     as.numeric(gsub("[^0-9]", "", as.character(x)))
##   }
##   x$Gross <- sapply(x$Gross, make_numeric)
##   x$TotalGross <- sapply(x$TotalGross, make_numeric)
##   x$Theaters <- sapply(x$Theaters, make_numeric)
## 
##   x$Change <- sapply(x$Change, function(x) as.numeric(gsub("[^-+0-9]", "", as.character(x))))
## 
##   ## add in weekend
##   x$weekend <- as.character(date)
##   x
## 
## }
## 
## x <- ymd("2011/01/07")
## all_weekends <- update(x, day=7 + 7*(0:51))
## l <- sapply(all_weekends, get_weekend, simplify=FALSE)
## names(l) <- unlist(Map(as.character, all_weekends))
## 
## movie_data_2011<- do.call(rbind,l)
## ##


## @knitr echo=FALSE
movie_data_2011 <- UsingR::movie_data_2011 ## remove when export issue is resolved


## @knitr 
x <- split(movie_data_2011, movie_data_2011$Movie)


## @knitr 
x <- Filter(function(d) is.na(d$Previous[1]), x)


## @knitr 

open_to_gross <- function(d) {
  ## return Genre and open to gross ratio
  list(movie = as.character(d$Movie[1]), 
       genre=as.character(d$Genre[1]), open=d$Gross[1], 
       total=max(d$TotalGross))
}
l <- sapply(x, open_to_gross, simplify=FALSE)


## @knitr 
d <- data.frame(l[[1]], stringsAsFactors=FALSE)
d <- Reduce(rbind, l[-1], d)


## @knitr 
d <- d[d$open > 5e5,]


## @knitr 
d$ratio <- d$total / d$open
summary(d$ratio)


## @knitr 
genres <- c("Comedy", "Action", "Adventure", "Horror")
d1 <- subset(d, subset=genre %in% genres)


## @knitr movie-ratio-by-genre, fig.keep="none"
boxplot(ratio ~ genre, data=d1) 


## @knitr echo=FALSE, out.width=singlewide
boxplot(ratio ~ genre, data=d1) 


## @knitr 
kruskal.test(ratio ~ factor(genre), d1)   


## @knitr 
names(morley)
oneway.test(Speed ~ Expt, data=morley)


## @knitr 
library(MASS)                           # load data set
oneway.test(MPG.highway ~ DriveTrain, data=Cars93)


## @knitr fig.keep="none"
boxplot(income ~ race, data=female.inc)      


## @knitr 
kruskal.test(income ~ race, data=female.inc)


## @knitr fig.keep="none"
boxplot(Driver.deaths ~ type, data=carsafety)


## @knitr 
oneway.test(Driver.deaths ~ type, data=carsafety,var.equal=FALSE)


## @knitr fig.keep="none"
boxplot(Other.deaths ~ type, data=carsafety)
oneway.test(Other.deaths ~ type, data=carsafety,var.equal=FALSE)


## @knitr fig.keep="none"
boxplot(BA ~ Hall.Fame.Membership, data=hall.fame)
oneway.test(BA ~ Hall.Fame.Membership, data=hall.fame, var.equal=FALSE)


## @knitr fig.keep="none"
lab1 <- c(4.13, 4.07, 4.04, 4.07, 4.05)
lab2 <- c(3.86, 3.85, 4.08, 4.11, 4.08)
lab3 <- c(4.00, 4.02, 4.01, 4.01, 4.04)
lab4 <- c(3.88, 3.89, 3.91, 3.96, 3.92)
chems <- stack(data.frame(lab1,lab2,lab3,lab4))
boxplot(values ~ ind, data=chems)       # var.equal unlikely
oneway.test(values ~ ind, data=chems, var.equal=FALSE)


## @knitr fig.keep="none"
type1 <- c(303, 293, 296, 299, 298)
type2 <- c(322, 326, 315, 318, 320, 320)
type3 <- c(309, 327, 317, 315)
wear <- list(type1=type1, type2= type2, type3=type3)
boxplot(wear)
oneway.test(values ~ ind, data=stack(wear))


## @knitr 
kruskal.test(weight ~ group, PlantGrowth)$p.value


## @knitr 
x <- c(63, 64, 95, 64, 60, 85)
y <- c(58, 56, 51, 84, 77)
z <- c(85, 79, 59, 89, 80, 71, 43)
d <- stack(list("test 1"=x, "test 2"=y, "test 3"=z))
kruskal.test(values ~ ind, data=d)


## @knitr 
oneway.test(values ~ ind, data=d)


## @knitr exam-boxplots, fig.keep="none"
x <- c(63, 64, 95, 64, 60, 85)
y <- c(58, 56, 51, 84, 77)
z <- c(85, 79, 59, 89, 80, 71, 43)
d <- stack(list("test 1"=x, "test 2"=y, "test 3"=z))
plot(values ~ ind, data=d, xlab = "test", ylab="grade")


## @knitr 
kruskal.test(values ~ ind, data=d)


## @knitr echo=FALSE, cache=FALSE
set.seed(3)


## @knitr cache=FALSE
mu1 <- 0; mu2 <- 1
x <- rnorm(15, mu1); y <- rnorm(15, mu2)
t.test(x, y, var.equal=TRUE)


## @knitr cache=FALSE
d <- stack(list(x=x, y=y))      # need named list.
res <- lm(values ~ ind, data=d)
short_summary(res)


## @knitr 
res <- lm(values ~ ind - 1, data = d)
short_summary(res)


## @knitr smoke-wt-boxplot, fig.keep="none"
library(UsingR)
d <- subset(babies, select=c("wt", "smoke"))   
plot(wt ~ factor(smoke), data=d,        # notice factor() for boxplot
     main="Birthweight by smoking level")


## @knitr echo=FALSE, out.width=singlewide, cache=FALSE
library(UsingR)
d <- subset(babies, select=c("wt", "smoke"))   
plot(wt ~ factor(smoke), data=d,        # notice factor() for boxplot
     main="Birthweight by smoking level")


## @knitr 
res <- lm(wt ~ factor(smoke), data=d)
short_summary(res)


## @knitr 
ewr.out <- subset(ewr, subset=inorout=="out", select=3:10)
out <- stack(ewr.out)
names(out) <- c("time","airline")
levels(out$airline)


## @knitr ewr-example, fig.keep="none"
plot(time ~ airline, data=out)
res <- lm(time ~ airline, data=out)
short_summary(res)


## @knitr ewr-example-2, fig.keep="none"
## trim airlines for display
out.trimmed <- subset(out, subset=airline %in% c("CO", "AA", "NW"))
res.aov.trimmed <- aov(time ~ airline, data=out.trimmed) 
TukeyHSD(res.aov.trimmed)
## 
res.aov <- aov(time ~ airline, data=out)
plot(TukeyHSD(res.aov), las=2)


## @knitr echo=FALSE, out.width=doublewide, cache=FALSE
plot(time ~ airline, data=out)
plot(TukeyHSD(res.aov), las=2)


## @knitr 
short_summary(lm(attendance ~ league, data=MLBattend))      


## @knitr 
short_summary(lm(y ~ limit, data = Traffic))
short_summary(lm(y ~ factor(year), data = Traffic))


## @knitr 
type1 <- c(303, 293, 296, 299, 298)
type2 <- c(322, 326, 315, 318, 320, 320)
type3 <- c(309, 327, 317, 315)
wear <- list(type1=type1, type2=type2, type3=type3)
wear.st <- stack(wear)


## @knitr 
oneway.test(values ~ ind, data = wear.st)$p.value


## @knitr 
summary(lm(values ~ ind, data = wear.st))


## @knitr 
oneway.test(values ~ ind, data = wear.st, var.equal=TRUE)$p.value


## @knitr 
short_summary(lm(mpg ~ factor(cyl), data=mtcars))


## @knitr 
short_summary(lm(mpg ~ factor(am), data=mtcars))    


## @knitr 
boxplot(amount ~ factor(year), data=npdb)


## @knitr 
boxplot(log(amount) ~ factor(year), data=npdb)      


## @knitr 
res <- lm(log(amount) ~ factor(year), subset= year < 2003, npdb)
short_summary(res)


## @knitr 
short_summary(lm(Speed ~ factor(Expt),data=morley))


## @knitr 
TukeyHSD(aov(Speed ~ factor(Expt), data=morley))


## @knitr 
res.aov <- aov(Other.deaths ~ type, data=carsafety)
TukeyHSD(res.aov)


## @knitr 
plot(count ~ spray, data=InsectSprays)
res.aov <- aov(count ~ spray, data=InsectSprays)
summary(res.aov)
plot(TukeyHSD(res.aov))  


## @knitr parallel-lines-smoker-take1, fig.keep="none"
babes <- subset(babies, subset = wt1 < 800)
plot(wt ~ wt1, data = babes, pch=as.numeric(smoke))


## @knitr 
res <- lm(wt ~ wt1 + factor(smoke), data = babes)
short_summary(res)


## @knitr parallel-lines-smoker, fig.keep="none"
plot(wt ~ wt1, pch=as.numeric(smoke), data=babes)
abline(107.0674, 0.1204)
abline(107.0674 - 8.3971, 0.1204, lty=2)


## @knitr echo=FALSE, out.width=singlewide
plot(wt ~ wt1, pch=as.numeric(smoke), data=babes)
abline(107.0674, 0.1204)
abline(107.0674 - 8.3971, 0.1204, lty=2)


## @knitr 
res.1 <- lm(wt ~ wt1, data = babes)
anova(res.1, res)


## @knitr cache=FALSE
res <- lm(time ~ age + gender, data=nym.2002)
short_summary(res)


## @knitr 
res = lm(mpg ~ wt + factor(am), data=mtcars)
short_summary(res)


## @knitr 
res <- lm(wt ~ gestation + wt1 + ht + factor(smoke), data=babies)
short_summary(res)


## @knitr 
kid.weights$BMI = with(kid.weights, (weight/2.2) / (height*2.54/100)^2)


## @knitr 
kid.weights$BMI = with(kid.weights, (weight/2.2) / (height*2.54/100)^2)
res.full = lm(BMI ~ age + gender, kid.weights)
res.age  = lm(BMI ~ age, kid.weights)
res.gender  = lm(BMI ~ gender, kid.weights)
anova(res.full, res.age)
anova(res.full, res.gender)


## @knitr 
require(MASS)                 
stepAIC(res)


## @knitr 
anova(lm(log(INCOME+1) ~ AGE + factor(EDUC), data=cfb))


## @knitr 
res = lm(temperature ~ hr + factor(gender), data=normtemp)
short_summary(res)


## @knitr interaction-plot-example, fig.keep="none", cache=FALSE
snacks$sugary <- cut(snacks$sugar, c(0, 20, 50, 100))
snacks$saturated <- cut(snacks$fat_sat, c(0, 7, 15, 25))
with(snacks, interaction.plot(saturated, sugary, response=calories))


## @knitr echo=FALSE, out.width=singlewide
snacks$sugary <- cut(snacks$sugar, c(0, 20, 50, 100))
snacks$saturated <- cut(snacks$fat_sat, c(0, 7, 15, 25))
with(snacks, interaction.plot(saturated, sugary, response=calories))


## @knitr echo=FALSE, results=FALSE
set.seed(50)
satisfied <- factor(sample(c("Yes", "No"), 12, TRUE))
tv <- factor(sample(c("0-1", "1-5","5+"), 12, TRUE))
x <- 0 + 12*(as.numeric(tv) - 1) + 6*(as.numeric(satisfied) -1)
y <- floor(pmax(0, rnorm(length(x), x, 1)))
y[y<=4] <- 0
d <- data.frame(x=y, satisfied=satisfied, tv = tv)


## @knitr 
d <- data.frame(
  x =c(0, 0, 16, 15, 14, 11, 26, 27, 23, 23, 36, 43),
  satisfied=c("yes", "no")[c(1,1,2,2,2,1,2,2,1,1,2,2)],
  tv=rep(c("0-1", "1-5", "5+"), c(5,3,4))
  )


## @knitr 
res.add <- lm(x ~ tv + satisfied, data=d)
res.int <- lm(x ~ tv * satisfied, data=d)


## @knitr 
summary(res.int)


## @knitr fig.keep="none"
with(d, interaction.plot( tv, satisfied, x))


## @knitr 
anova(res.int, res.add)


## @knitr 
short_summary(res.add)


## @knitr cache=FALSE
res.int <- lm(calories ~ saturated * sugary, snacks)
res.add <- lm(calories ~ saturated + sugary, snacks)
anova(res.int, res.add)


## @knitr 
library(MASS)
out <- t.test(shoes$A, shoes$B, paired=TRUE)
out$p.value


## @knitr 
d <- stack(shoes)
names(d) <- c("value", "sole_type")
d$person <- gl(10, 1, 20, labels=LETTERS[1:10])
xtabs(value ~ sole_type + person, data=d)

res <- aov(value ~ sole_type + person, data=d)
summary(res)


## @knitr 
summary(aov(value ~ sole_type, data=d))


## @knitr 
d <- data.frame(
  percent_change = c(1.3, -0.7,  0.3, -1.4, 0.5 ,0.9 ,1.6, 1.2, 
    0.3, 1.8, 1.7, 1.3),
  ethnicity = gl(4, 1, 12, labels=paste("Ethic", 1:4)),
  trt = gl(3, 4, 12, labels=c("Ctrl", "Soda", "Vitamin"))
  )
xtabs(percent_change ~ trt + ethnicity, data=d)


## @knitr 
res <- aov(percent_change ~ trt + ethnicity, data=d)
coef(res)[1:3]                          # only intercept, trt


## @knitr 
summary(res)


## @knitr 
likability <- c(-1, 4, 0, -1, 4, 1, 6, 2, 7, 1, 2, 2, 7, 5, 2, 3, 6, 1)
web <- gl(2, 9, 18, c("N", "Y"))
tv <- gl(3, 6, 18, labels=c(0, "1-2", "3+"))

res.web  <- lm(likability ~ web)
res.tv   <- lm(likability ~ tv)
res.both <- lm(likability ~ tv + web)


## @knitr 
summary(res.web)


## @knitr 
anova(res.tv, res.both)


## @knitr eval=FALSE
## ftable(xtabs(UBP ~ person + replicate + grip.type, data=grip))


## @knitr 
with(grip, interaction.plot(grip.type,person,UBP))
res.int   <- lm(UBP ~ person * grip.type, data=grip)
res.noint <- lm(UBP ~ person + grip.type, data=grip)
res.per   <- lm(UBP ~ person,             data=grip)
res.grip  <- lm(UBP ~ grip.type,          data=grip)
res.none  <- lm(UBP ~ 1,                  data=grip)


## @knitr 
anova(res.noint,res.int)


## @knitr 
anova(res.none,res.per)


## @knitr 
anova(res.none,res.grip)


## @knitr 
library(MASS, verbose=FALSE)
stepAIC(res.int) 


## @knitr 
res.full = lm(mpg ~ factor(cyl) * factor(am), data=mtcars)
short_summary(res.full)


## @knitr 
res.add <- lm(mpg ~ factor(cyl) + factor(am), data=mtcars)
res.cyl <- lm(mpg ~ factor(cyl),              data=mtcars)
res.am <- lm(mpg ~ factor(am),                data=mtcars)
anova(res.am, res.add)
anova(res.cyl, res.add)


## @knitr 
with(ToothGrowth, interaction.plot(dose, supp, len))
res.full <- lm(len ~ supp + dose, data=ToothGrowth)
res.add  <- lm(len ~ supp * dose, data=ToothGrowth)
anova(res.full, res.add)


## @knitr fig.keep="none"
with(OrchardSprays, interaction.plot(rowpos,treatment,decrease))


## @knitr 
res.full <- lm(decrease ~ rowpos * treatment, data=OrchardSprays)
res.add  <- lm(decrease ~ rowpos + treatment, data=OrchardSprays)
anova(res.full, res.add)


## @knitr fig.keep="none"
d <- data.frame(
  rating=c(8, 6, 8, 4),
  food=gl(2, 2, 4, labels=c("bread", "corn")),
  butter=gl(2, 1, 4, labels=c("yes", "no")))
xtabs(rating ~ butter + food, d)
with(d, interaction.plot(butter, food, rating)) # not shown


## @knitr eval=FALSE
## summary(lm(rating ~ butter * food, d))


