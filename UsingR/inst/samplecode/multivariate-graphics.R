
mapply <- function(...) invisible(base:::mapply(...))


## @knitr iris-scatter, eval=FALSE
## with(iris,
##      plot(Sepal.Length, Sepal.Width,
##           pch=as.numeric(Species), cex=1.2))
## legend(6.1, 4.4, c("setosa", "versicolor", "virginica"),
##        cex=1.5, pch=1:3)


## @knitr iris-scatter-abline, eval=FALSE
## fm <- Sepal.Width ~ Sepal.Length
## plot(fm, iris, pch=as.numeric(Species))
## out <- mapply(function(i, x) abline(lm(fm, data=x), lty=i),
##               i=1:3, x = split(iris, iris$Species))
## legend(6.5, 4.4, levels(iris$Species), cex=1.5, lty=1:3)


## @knitr echo=FALSE, out.width=doublewide
with(iris,
     plot(Sepal.Length, Sepal.Width,
          pch=as.numeric(Species), cex=1.2))
legend(6.1, 4.4, c("setosa", "versicolor", "virginica"),
       cex=1.5, pch=1:3)
fm <- Sepal.Width ~ Sepal.Length
plot(fm, iris, pch=as.numeric(Species))
out <- mapply(function(i, x) abline(lm(fm, data=x), lty=i),
              i=1:3, x = split(iris, iris$Species))
legend(6.5, 4.4, levels(iris$Species), cex=1.5, lty=1:3)


## @knitr cars93-highway-weight-price, eval=FALSE
## Cars93 <- transform(Cars93, price=cut(Price, c(0, 15, 30, 75),
##                    labels=c("cheap", "affordable", "expensive")))
## plot(MPG.highway ~ Weight, Cars93, pch=as.numeric(price))
## legend(3500, 50, levels(Cars93$price), pch=1:3)


## @knitr echo=FALSE, out.width=singlewide
Cars93 <- transform(Cars93, price=cut(Price, c(0, 15, 30, 75), 
                   labels=c("cheap", "affordable", "expensive")))
plot(MPG.highway ~ Weight, Cars93, pch=as.numeric(price))
legend(3500, 50, levels(Cars93$price), pch=1:3)


## @knitr cars93-highway-weight-price-par, eval=FALSE
## l <- split(Cars93, Cars93$price)
## ## 3 graphics in one
## par(mfrow=c(1,length(l)))
## ## common to each graphic
## fm <- MPG.highway ~ Weight
## xlim <- range(Cars93$Weight)
## ylim <- range(Cars93$MPG.highway)
## 
## mapply(function(x, nm) {
##   plot(fm, data=x, main=nm, xlim=xlim, ylim=ylim)
##   abline(lm(fm, data=x))
## }, l, names(l))


## @knitr echo=FALSE, out.width=triplewide
l <- split(Cars93, Cars93$price)
## 3 graphics in one
##par(mfrow=c(1,3))
## common to each graphic
f <- MPG.highway ~ Weight
xlim <- range(Cars93$Weight)
ylim <- range(Cars93$MPG.highway)

mapply(function(x, nm) {
  plot(f, x, main=nm, xlim=xlim, ylim=ylim)
  abline(lm(f, x))
}, l, names(l))


## @knitr SAT-scores-by-salary-percent, eval=FALSE
## plot(total ~ salary, data=SAT, cex=sqrt(perc/10),
##      pch=16,                                       # filled circles
##      col=rgb(red=0, green=0, blue=0, alpha=0.250)) # use alpha


## @knitr echo=FALSE, out.width=singlewide
plot(total ~ salary, data=SAT, cex=sqrt(perc/10), 
     pch=16,                                       # filled circles
     col=rgb(red=0, green=0, blue=0, alpha=0.250)) # use alpha


## @knitr pairs-plot-iris, eval=FALSE
## species <- iris$Species
## values <- Filter(is.numeric, iris)
## pairs(values, col=species)


## @knitr echo=FALSE, out.width=singlewide
species <- iris$Species
values <- Filter(is.numeric, iris)
pairs(values, col=species)


## @knitr parallel-coord, eval=FALSE
## x <- state.x77
## case <- "New York"
## ind <- rownames(x) == case
## parcoord(state.x77,
##          col=gray(c(.8, .2))[ind + 1],  # light/dark gray
##          lwd=(1:2)[ind + 1], las = 2)


## @knitr echo=FALSE, fig.width=7, fig.height=3.5
x <- state.x77
case <- "New York"
ind <- rownames(x) == case
parcoord(state.x77, 
         col=gray(c(.8, .2))[ind + 1],  # light/dark gray
         lwd=(1:2)[ind + 1], las = 2)


## @knitr 
x <- sapply(as.data.frame(state.x77), rank)
rownames(x) <- rownames(state.x77)


## @knitr heatmap, eval=FALSE
## heatmap(x, Rowv=NA, Colv=NA,
##         scale="column",                 # scale columns
##         margins=c(8,6),                 # leave room for labels
##         col=rev(gray.colors(50)))       # darker -> larger
## 


## @knitr echo=FALSE, out.width=singlewide
heatmap(x, Rowv=NA, Colv=NA, 
        scale="column",                 # scale columns
        margins=c(8,6),                 # leave room for labels
        col=rev(gray.colors(50)))       # darker -> larger



## @knitr fig.keep="none"
plot(calories ~ sugars, pch=shelf, data=UScereal)
legend(0, 400, legend=paste("shelf", 1:3),  pch=1:3)


## @knitr 
values <- Filter(is.numeric, iris)
cols <- c("red", "blue", "green")
cols <- cols[as.numeric(iris$Species)]
parcoord(values, col=cols)


## @knitr 
require(MASS)
d <- Filter(is.numeric, UScereal)
pairs(d)


## @knitr 
with(UScereal, c(calories=cor(fat, calories), sugars=cor(fat, sugars)))


## @knitr fig.keep="none"
 plot(HR ~ H, batting, cex=sqrt(SO)/10)


## @knitr eval=FALSE
## xyplot(MPG.highway ~ Weight | Price, data=Cars93)


## @knitr eval=FALSE
## Cars93 = transform(Cars93, price=cut(Price, c(0, 15, 30, 75),
##                    labels=c("cheap", "affordable", "expensive")))
## xyplot(MPG.highway ~ Weight | price, data=Cars93)


## @knitr 
prices <- equal.count(Cars93$Price, number=3, overlap=0)


## @knitr 
panel=function(x, y) {
  panel.xyplot(x, y)
  panel.lmline(x, y)
}


## @knitr xyplot-highway-weight-price, eval=FALSE
## xyplot(MPG.highway ~ Weight | price, data=Cars93,
##        layout=c(3,1),
##        type=c("p", "r")
## )


## @knitr echo=FALSE, fig.height=4.5, out.width=singlewide, cache=FALSE
xyplot(MPG.highway ~ Weight | price, data=Cars93,
       layout=c(3,1),
       type=c("p", "r")
)


## @knitr dotplot-smoke-wt, eval=FALSE
## dotplot(factor(smoke) ~ wt, data=babies, subset=wt < 999 & ded==3)


## @knitr bwplot-smoke-wt, eval=FALSE
## bwplot(factor(smoke) ~ wt, data=babies, subset=wt < 999)


## @knitr dotplot-smoke-wt-ded, eval=FALSE
## dotplot(factor(smoke) ~ wt | factor(ded), data=babies, subset=wt < 999)


## @knitr histogram-wt-smoke, eval=FALSE
## histogram(~ wt | factor(smoke), data=babies, subset=wt < 999)


## @knitr densityplot-wt-smoke, eval=FALSE
## densityplot( ~ wt | factor(smoke), data=babies, subset=wt < 999)


## @knitr fig.keep="none"
require(MASS)
require(lattice)
dotplot(Expt ~ Speed, data=michelson)


## @knitr fig.keep="none"
require(lattice)
bwplot(teamID ~ I(H/AB), batting)


## @knitr 
sapply(split(batting, batting$teamID), 
       function(DF) with(DF, median(H/AB)))


## @knitr 
p <- ggplot(Cars93)


## @knitr ggplot-weight-mpg-heighway
p <- p + aes(x=Weight, y = MPG.highway, color=Origin)


## @knitr aes
aes(x=sqrt(Weight/1000), y=MPG.highway, color=Origin)


## @knitr ggplot-weight-highway, echo=FALSE, eval=FALSE
## p <- ggplot(Cars93, aes(x=sqrt(Weight/1000), y=MPG.highway,
##                         color=Origin)) + geom_point(cex=3)
## p                                       # show the plot


## @knitr ggplot-two-geoms, echo=FALSE, eval=FALSE
## f <- function(x) x^2
## x <- -2:2; y <- f(x)
## p <- ggplot(data.frame(x=x, y=y), aes(x=x, y=y)) +
##   geom_line(color="black", alpha=0.25) +      # set color attribute
##   geom_point(pch=16, cex=5)           # set plot character and size
## p


## @knitr ggplot2-basics, echo=FALSE, out.width=doublewide
p <- ggplot(Cars93, aes(x=sqrt(Weight/1000), y=MPG.highway, 
                        color=Origin)) + geom_point(cex=3)
p                                       # show the plot
f <- function(x) x^2
x <- -2:2; y <- f(x)
p <- ggplot(data.frame(x=x, y=y), aes(x=x, y=y)) + 
  geom_line(color="black", alpha=0.25) +      # set color attribute
  geom_point(pch=16, cex=5)           # set plot character and size
p


## @knitr eval=FALSE
## p <- ggplot(Cars93, aes(x=sqrt(Weight/1000), y=MPG.highway, 
##                         color=Origin)) + geom_point(cex=3)
## p                                       # show the plot


## @knitr ggplot-curve-like, fig.keep="none"
f <- function(x) x^2
x <- seq(-2, 2, length=100)
p <- ggplot(data.frame(x=x, y=f(x)), aes(x=x, y=y))
p + geom_line()                         # not shown


## @knitr eval=FALSE
## f <- function(x) x^2
## x <- -2:2; y <- f(x)
## p <- ggplot(data.frame(x=x, y=y), aes(x=x, y=y)) + 
##   geom_line(color="black", alpha=0.25) +      # set color attribute
##   geom_point(pch=16, cex=5)           # set plot character and size
## p


## @knitr ggplot-boxplot-cylinder-highway, echo=FALSE, eval=FALSE
## p <- ggplot(Cars93, aes(x=Cylinders, y=MPG.highway))
## p + geom_boxplot()


## @knitr ggplot-boxplot-cylinder-highway-jitter, echo=FALSE, eval=FALSE
## p + geom_boxplot() +
##   geom_jitter(position=position_jitter(w = 0.1), alpha=.25)


## @knitr echo=FALSE, out.width=doublewide
p <- ggplot(Cars93, aes(x=Cylinders, y=MPG.highway))
p + geom_boxplot()
p + geom_boxplot() + 
  geom_jitter(position=position_jitter(w = 0.1), alpha=.25)


## @knitr eval=FALSE
## p <- ggplot(Cars93, aes(x=Cylinders, y=MPG.highway))
## p + geom_boxplot()


## @knitr eval=FALSE
## p + geom_boxplot() + 
##   geom_jitter(position=position_jitter(w = 0.1), alpha=.25)


## @knitr morley-speed-run-by-expt, eval=FALSE
## m <- subset(morley, Expt %in% 1:2)      # just first two experiments
## p <- ggplot(m, aes(x=Run, y=Speed, group=Expt, color=Expt))
## p + geom_line()                         # connect x and y with lines


## @knitr echo=FALSE, out.width=singlewide
m <- subset(morley, Expt %in% 1:2)      # just first two experiments
p <- ggplot(m, aes(x=Run, y=Speed, group=Expt, color=Expt))
p + geom_line()                         # connect x and y with lines


## @knitr ggplot-hist-highway, eval=FALSE
## p <- ggplot(Cars93, aes(x = MPG.highway))
## p + stat_bin(binwidth=5)


## @knitr ggplot-hist-highway-density, eval=FALSE
## p <- ggplot(Cars93, aes(x = MPG.highway, y=..density..)) +
##   stat_bin(binwidth=5)
## p


## @knitr echo=FALSE, out.width=doublewide
p <- ggplot(Cars93, aes(x = MPG.highway))
p + stat_bin(binwidth=5)
p <- ggplot(Cars93, aes(x = MPG.highway, y=..density..)) + 
  stat_bin(binwidth=5)
p


## @knitr eval=FALSE
## p <- ggplot(Cars93, aes(x=MPG.highway, y=..density..)) # scale y
## p + geom_histogram(alpha=0.5) + geom_density()


## @knitr eval=FALSE
## p <- ggplot(Cars93, aes(x=Weight, y=MPG.highway)) + geom_point()
## p + geom_smooth()


## @knitr ggplot-stat-smooth, echo=FALSE, eval=FALSE
## ## hide warning for method="loess"
## p <- ggplot(Cars93, aes(x=Weight, y=MPG.highway)) + geom_point()
## p + geom_smooth(method="loess")


## @knitr ggplot-stat-smooth-lm, eval=FALSE
## p + geom_smooth(method="lm", se=FALSE)


## @knitr ggplot-stat-smooth-lm-poly, eval=FALSE
## p + geom_smooth(method="lm", formula=y ~ poly(x,2), se=FALSE)


## @knitr echo=FALSE, out.width=triplewide
## hide warning for method="loess"
p <- ggplot(Cars93, aes(x=Weight, y=MPG.highway)) + geom_point()
p + geom_smooth(method="loess")
p + geom_smooth(method="lm", se=FALSE)
p + geom_smooth(method="lm", formula=y ~ poly(x,2), se=FALSE)


## @knitr highway-weight-ggplot, echo=FALSE,  fig.height=4.5, out.width=singlewide
Cars93 <- transform(Cars93, price=cut(Price, c(0, 15, 30, 75), 
                   labels=c("cheap", "affordable", "expensive")))
p <- ggplot(Cars93) + aes(x=Weight, y = MPG.highway) + 
  geom_point(cex=3) + 
  geom_smooth(method="lm", se=FALSE)  +
  facet_grid( ~ price)
p


## @knitr ggplot-facet-wrap, eval=FALSE
## p <- ggplot(PearsonLee, aes(y=child,x=parent))
## p + geom_point(alpha=0.5) + geom_smooth(method="loess") +
##   facet_grid(par ~ chl)


## @knitr eval=FALSE
## Cars93 <- transform(Cars93, price=cut(Price, c(0, 15, 30, 75), 
##                    labels=c("cheap", "affordable", "expensive")))
## p <- ggplot(Cars93) + aes(x=Weight, y = MPG.highway) + 
##   geom_point(cex=3) + 
##   geom_smooth(method="lm", se=FALSE)  +
##   facet_grid( ~ price)
## p


## @knitr ggplot-facet-wrap-grid, eval=FALSE
## p <-  ggplot(PearsonLee, aes(y=child,x=parent))
## p + geom_point(alpha=0.5) +
##   geom_smooth(method="loess") +
##   facet_grid(par ~ chl, margins="chl")  # or margins=TRUE for both


## @knitr echo=FALSE, out.width=doublewide, warning=FALSE
p <- ggplot(PearsonLee, aes(y=child,x=parent))
p + geom_point(alpha=0.5) + geom_smooth(method="loess") + 
  facet_grid(par ~ chl)
p <-  ggplot(PearsonLee, aes(y=child,x=parent))
p + geom_point(alpha=0.5) + 
  geom_smooth(method="loess") + 
  facet_grid(par ~ chl, margins="chl")  # or margins=TRUE for both 


## @knitr ggplot-hist-facet-wrap, eval=FALSE
## p <- ggplot(morley, aes(x=Speed)) + geom_histogram(binwidth=50)
## p + facet_wrap( ~ Expt)


## @knitr echo=FALSE, out.width=singlewide, warning=FALSE
p <- ggplot(morley, aes(x=Speed)) + geom_histogram(binwidth=50)
p + facet_wrap( ~ Expt)


## @knitr 
p <- ggplot(mtcars)
p + aes(x=factor(gear), y=mpg) + geom_boxplot() # boxplots
p + aes(x=wt, y=mpg) + geom_point() + geom_smooth(method="lm") # scatter plot
p + aes(x=hp, y=mpg, col=factor(am), pch=factor(am)) + # add color and shape for transmission
  geom_point() + 
  facet_grid(gear ~ cyl)


