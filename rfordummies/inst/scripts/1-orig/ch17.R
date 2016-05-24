# Chapter 17 - Creating Faceted Graphics with Lattice


# Creating a Lattice Plot

str(mtcars)

## Loading the lattice package

library("lattice")

## Making a lattice scatterplot

xyplot(mpg ~ hp | factor(cyl), data=mtcars)

## Adding trend lines

xyplot(mpg ~ hp | factor(cyl), data=mtcars,
   type=c("p", "r"))


# Changing Plot Options

## Adding titles and labels

xyplot(mpg ~ hp | factor(cyl), data=mtcars,
   type=c("p", "r"),
   main="Fuel economy vs. Performance",
   xlab="Performance (horse power)",
   ylab="Fuel economy (miles per gallon)",
)

xyplot(mpg ~ hp | factor(cyl), data=mtcars,
   type=c("p", "r"),
   main=list(
       label="Fuel economy vs. Performance given Number of Cylinders",
       cex=0.75)
)

## Changing the font size of titles and labels

xyplot(mpg ~ hp | factor(cyl), data=mtcars,
   type=c("p", "r"),
   main=list(
       label="Fuel economy vs. Performance given Number of Cylinders",
       cex=0.75),
   xlab=list(
       label="Performance (horse power)",
       cex=0.75),
   ylab=list(
       label="Fuel economy (miles per gallon)",
       cex=0.75),
   scales=list(cex=0.5)
)


## Using themes to modify plot options

xyplot(mpg ~ hp | factor(cyl), data=mtcars,
   type=c("p", "r"),
   par.settings=simpleTheme(col="red", col.line="blue")
)


# Plotting Different Types

## Making a bar chart

mtcars$cars <- rownames(mtcars)

barchart(cars ~ mpg | factor(cyl), data=mtcars,
   main="barchart",
   scales=list(cex=0.5),
   layout=c(3, 1)
)

## Making a box-and-whisker plot

bwplot(~ hp | factor(cyl), data=mtcars, main="bwplot")


# Plotting Data in Groups

## Using data in tall format

str(longley)
library("reshape2")
mlongley <- melt(longley, id.vars="Year")

str(mlongley)

xyplot(value ~ Year | variable, data=mlongley,
   layout=c(6, 1),
   par.strip.text=list(cex=0.7),
   scales=list(cex=0.7)
)

## Creating a chart with groups

mtcars$cars <- rownames(mtcars)
mtcars$am <- with(mtcars, ifelse(am==0, "Automatic", "Manual"))

barchart(cars ~ mpg | factor(cyl), data=mtcars,
   group=am,
   scales=list(cex=0.5),
   layout=c(3, 1),
)

## Adding a key

barchart(cars ~ mpg | factor(cyl), data=mtcars,
   main="barchart with groups",
   group=am,
   auto.key=TRUE,
   par.settings = simpleTheme(col=c("grey80", "grey20")),
   scales=list(cex=0.5),
   layout=c(3, 1)
)


# Printing and Saving a Lattice Plot

## Assigning a lattice plot to an object

my.plot <- xyplot(mpg ~ hp | cyl, data=mtcars)
class(my.plot)

## Printing a lattice plot in a script

xyplot(mpg ~ hp | cyl, data=mtcars)

my.plot <- xyplot(mpg ~ hp | cyl, data=mtcars)
print(my.plot)


## Saving a lattice plot to file

setwd("~/")
trellis.device(device="png", filename="xyplot.png")
print(my.plot)
dev.off()

