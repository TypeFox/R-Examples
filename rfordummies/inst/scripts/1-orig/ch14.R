# Chapter 14
# Summarizing Data

# Starting with the Right Data

## Using factors or numeric data

## Counting unique values

sapply(mtcars, function(x) length(unique(x)))

## Preparing the data

cars <- mtcars[c(1,2,9,10)]
cars$gear <- ordered(cars$gear)
cars$am <- factor(cars$am, labels=c('auto', 'manual'))
str(cars)

# Describing Continuous Variables

## Talking about the center of your data

mean(cars$mpg)
median(cars$cyl)

## Describing the variation
sd(cars$mpg)

## Checking the quantiles

### Calculating the range
range(cars$mpg)

### Calculating the quantiles
quantile(cars$mpg)

### Getting on speed with the quantile function
quantile(cars$mpg, probs=c(0.05, 0.95))

# Describing Categories

## Counting appearances

### Creating a table
amtable <- table(cars$am)
amtable

### Working with tables

## Calculating proportions
amtable/sum(amtable)
prop.table(amtable)

## Finding the center
id <- amtable == max(amtable)
names(amtable)[id]

# Describing Distributions

## Plotting histograms

### Making the plot
hist(cars$mpg, col='grey')

### Playing with breaks
hist(cars$mpg, breaks=c(5,15,25,35))

## Using frequencies or densities

### Creating a density plot

mpgdens <- density(cars$mpg)
plot(mpgdens)

### Plotting densities in a histogram
hist(cars$mpg, col='grey', freq=FALSE)
lines(mpgdens)

# Describing Multiple Variables

## Summarizing a complete dataset

### Getting the output
summary(cars)

### Fixing a problem

cars$cyl <- as.factor(cars$cyl)

## Plotting quantiles for subgroups

boxplot(mpg ~ cyl, data=cars)

## Tracking correlations

names(iris)

### Looking at relations
plot(iris[-5])

### Getting the numbers

with(iris, cor(Petal.Width, Petal.Length))

### Calculating correlations for multiple variables

iris.cor <- cor(iris[-5])
str(iris.cor)

iris.cor['Petal.Width', 'Petal.Length']

### Dealing with missing values

# Working with Tables

## Creating a two-way table

### Creating a table from two variables

with(cars, table(am, gear))

### Creating tables from a matrix

trial <- matrix(c(34,11,9,32), ncol=2)
colnames(trial) <- c('sick', 'healthy')
rownames(trial) <- c('risk', 'no_risk')
trial.table <- as.table(trial)
trial.table

### Extracting the numbers

trial.table['risk', 'sick']

##Converting tables to a data frame

trial.df <- as.data.frame(trial)
str(trial.df)

trial.table.df <- as.data.frame(trial.table)
str(trial.table.df)

## Looking at margins and proportions

### Adding margins to the table

addmargins(trial.table)
addmargins(trial.table,margin=2)

### Calculating proportions

prop.table(trial.table)

### Calculating proportions over columns and rows
prop.table(trial.table, margin=1)
