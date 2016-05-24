# Chapter 13 - Manipulating and Processing Data

# Deciding on the Most Appropriate Data Structure

# Creating Subsets of Your Data

## Understanding the three subset operators
## Understanding the five ways of specifying the subset

str(islands)
islands[]
islands[c(8, 1, 1, 42)]
islands[-(3:46)]
islands[islands < 20]
islands[c("Madagascar", "Cuba")]

## Subsetting data frames

str(iris)
iris[1:5, ]
iris[, c("Sepal.Length", "Sepal.Width")]
iris[, 'Sepal.Length']
iris[, 'Sepal.Length', drop=FALSE]
iris['Sepal.Length']
iris[1:5, c("Sepal.Length", "Sepal.Width")]

### Taking samples from data

sample(1:6, 10, replace=TRUE)

set.seed(1)
sample(1:6, 10, replace=TRUE)
sample(1:6, 10, replace=TRUE)

set.seed(1)
sample(1:6, 10, replace=TRUE)

set.seed(123)
index <- sample(1:nrow(iris), 5)
index
iris[index, ]

### Removing duplicate data

duplicated(c(1,2,1,3,1,4))
duplicated(iris)
which(duplicated(iris))
iris[!duplicated(iris), ]

index <- which(duplicated(iris))
iris[-index, ]

### Removing rows with missing data

str(airquality)
complete.cases(airquality)

x <- airquality[complete.cases(airquality), ]
str(x)
x <- na.omit(airquality)



# Adding Calculated Fields to Data

## Doing arithmetic on columns of a data frame

x <- iris$Sepal.Length / iris$Sepal.Width
head(x)

## Using with and within to improve code readability

y <- with(iris, Sepal.Length / Sepal.Width)
head(y)
identical(x, y)

iris$ratio <- iris$Sepal.Length / iris$Sepal.Width
iris <- within(iris, ratio <- Sepal.Length / Sepal.Width)
head(iris$ratio)

## Creating subgroups or bins of data

### Using cut to create a fixed number of subgroups

head(state.x77)
frost <- state.x77[, "Frost"]
head(frost, 5)
cut(frost, 3, include.lowest=TRUE)

### Adding labels to cut

cut(frost, 3, include.lowest=TRUE, labels=c("Low", "Med", "High"))

### Using table to count the number of observations

x <- cut(frost, 3, include.lowest=TRUE, labels=c("Low", "Med", "High"))
table(x)
x


# Combining and Merging Data Sets

## Creating sample data to illustrate merging

all.states <- as.data.frame(state.x77)
all.states$Name <- rownames(state.x77)
rownames(all.states) <- NULL
str(all.states)

### Creating a subset of cold states

cold.states <- all.states[all.states$Frost>150, c("Name", "Frost")]
cold.states

### Creating a subset of large states

large.states <- all.states[all.states$Area>=100000, c("Name", "Area")]
large.states

## Using the merge() function

### Using merge to find the intersection of data

merge(cold.states, large.states)

### Understanding the different types of merge

merge(cold.states, large.states, all=TRUE)


## Working with lookup tables

### Finding a match

index <- match(cold.states$Name, large.states$Name)
index

large.states[na.omit(index), ]

### Making sense of %in%

index <- cold.states$Name %in% large.states$Name
index
!is.na(match(cold.states$Name,large.states$Name))
cold.states[index, ]

# Sorting and Ordering Data

some.states <- data.frame(
     Region = state.region,
     state.x77)

some.states <- some.states[1:10, 1:3]
some.states

## Sorting vectors

### Sorting a vector in ascending order

sort(some.states$Population)

### Sorting a vector in decreasing order

sort(some.states$Population, decreasing=TRUE)

## Sorting data frames

### Getting the order

order.pop <- order(some.states$Population)
order.pop

some.states$Population[order.pop]

## Sorting a data frame in ascending order

some.states[order.pop, ]
order(some.states$Population)
order(some.states$Population, decreasing=TRUE)

some.states[order(some.states$Population, decreasing=TRUE), ]

### Sorting on more than one column

index <- with(some.states, order(Region, Population))
some.states[index, ]

### Sorting multiple columns in mixed order
index <- order(-xtfrm(some.states$Region), some.states$Population)
some.states[index, ]

# Traversing Your Data with the Apply Functions

## Using the apply() function to summarize arrays

str(Titanic)
apply(Titanic, 1, sum)
apply(Titanic, 3, sum)
apply(Titanic, c(3, 4), sum)

## Using lapply() and sapply() to traverse a list or data frame

lapply(iris, class)
sapply(iris, class)
sapply(iris, mean)
sapply(iris, function(x) ifelse(is.numeric(x), mean(x), NA))

## Using tapply() to create tabular summaries

tapply(iris$Sepal.Length, iris$Species, mean)
with(iris, tapply(Sepal.Length, Species, mean))

### Using tapply() to create higher-dimensional tables

str(mtcars)
cars <- within(mtcars,
    am <- factor(am, levels=0:1, labels=c("Automatic", "Manual"))
)

with(cars, tapply(mpg, am, mean))
with(cars, tapply(mpg, list(gear, am), mean))

### Using aggregate()

with(cars, aggregate(mpg, list(gear=gear, am=am), mean))

# Getting to Know the Formula Interface


aggregate(mpg ~ gear + am, data=cars, mean)

aov(mpg ~ gear + am, data=cars)

library(lattice)
xyplot(mpg ~ gear + am, data=cars)


# Whipping Your Data into Shape


## Understanding data in long and wide format


## Getting started with the reshape2 package

install.packages("reshape2")
library("reshape2")

goals <- data.frame(
    Game = c("1st", "2nd", "3rd", "4th"),
    Venue = c("Bruges", "Ghent", "Ghent", "Bruges"),
    Granny = c(12, 4, 5, 6),
    Geraldine = c(5, 4, 2, 4),
    Gertrude = c(11, 5, 6, 7)
)

## Melting data to long format

mgoals <- melt(goals)
mgoals <- melt(goals, id.vars=c("Game", "Venue"))
mgoals

## Casting data to wide format

dcast(mgoals,  Venue + Game ~ variable, sum)
dcast(mgoals, variable ~ Venue , sum)
dcast(mgoals,  Venue ~ variable , sum)

dcast(mgoals,  Venue + variable ~ Game , sum)

library(ggplot2)
ggplot(mgoals, aes(x=variable, y=value, fill=Game)) + geom_bar(stat="identity")
