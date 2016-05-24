## ----cars_setup----------------------------------------------------------
# ggplot2 examples
library(ggplot2)
#use color brewer as default discrete colors
scale_colour_discrete <- function(...) scale_color_brewer(palette="Set1", ...)
scale_fill_discrete <- function(...) scale_fill_brewer(palette="Set1", ...)

data('mtcars')
# create factors with value labels
mtcars$gear <- factor(mtcars$gear,levels=c(3,4,5),
     labels=c("3gears","4gears","5gears"))
mtcars$am <- factor(mtcars$am,levels=c(0,1),
     labels=c("Automatic","Manual"))
mtcars$cyl <- factor(mtcars$cyl,levels=c(4,6,8),
   labels=c("4cyl","6cyl","8cyl"))
head(mtcars)

## ----cars_density, dev='png', warning=FALSE------------------------------
qplot(mpg, data=mtcars, geom="density", fill=gear, alpha=I(.5),
   main="Distribution of Gas Milage", xlab="Miles Per Gallon",
   ylab="Density")

## ----cars_scatter, dev='png', warning=FALSE------------------------------
qplot(hp, mpg, data=mtcars, shape=am, color=am,
   facets=gear~cyl, size=I(3),
   xlab="Horsepower", ylab="Miles per Gallon")

## ----cars_regressions, dev='png', warning=FALSE--------------------------
ggplot(data = mtcars, aes(x = wt, y = mpg, color = cyl)) + geom_point() +
  geom_smooth(method="lm") +
  labs(main="Regression of MPG on Weight",
   xlab="Weight", ylab="Miles per Gallon")

## ----cars_boxplots, dev='png', warning=FALSE-----------------------------
qplot(gear, mpg, data=mtcars, geom=c("boxplot", "jitter"),
   fill=gear, main="Mileage by Gear Number",
   xlab="", ylab="Miles per Gallon")

