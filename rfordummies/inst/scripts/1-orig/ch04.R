# Chapter 4 - Getting Started with Arithmetic

# Working with Numbers, Infinity, and Missing Values

## Doing basic arithmetic

### Using arithmetic operators


baskets.of.Granny <- c(12,4,4,6,9,3)
baskets.of.Geraldine <- c(5,3,2,2,12,9)


Granny.money <- baskets.of.Granny * 120
Geraldine.money <- baskets.of.Geraldine * 145


Granny.money + Geraldine.money


baskets.of.Granny * 120 + baskets.of.Geraldine * 145

### Controlling the order of the operations
4 + 2 * 3
(4 + 2)* 3

## Using mathematical functions

### Calculating logarithms and exponentials

log(1:3)
log(1:3,base=6)

x <- log(1:3)
exp(x)

### Putting the science in scientific notation
1.33e4

4.12e-2

1.2e6 / 2e3


### Rounding numbers

round(123.456,digits=2)
round(-123.456,digits=-2)
signif(-123.456,digits=4)

### Using trigonometric functions

cos(120)
cos(120*pi/180)

## Calculating whole vectors

`+`(2,3)


##To infinity and beyond

### Using infinity

2/0
4 - Inf
is.finite(10^(305:310))

### Dealing with undefined outcomes
Inf / Inf
NaN + 4


### Dealing with missing values

x <- NA
x + 4

log(x)

is.na(x)

### Calculating infinite, undefined, and missing values


# Organizing Data in Vectors

## Discovering the properties of vectors

### Looking at the structure of a vector

str(baskets.of.Granny)
length(baskets.of.Granny)
authors <- c("Andrie", "Joris")
str(authors)

### Testing vector types




is.numeric(baskets.of.Granny)
is.integer(baskets.of.Granny)

x <- c(4L,6L)
is.integer(x)

## Creating vectors

seq(from = 4.5, to = 2.5, by = -0.5)


seq(from = -2.7, to = 1.3, length.out = 9)


baskets.of.Granny <- c(12,4,4,6,9,3)
baskets.of.Geraldine <- c(5,3,2,2,12,9)

## Combining vectors

all.baskets <-c(baskets.of.Granny, baskets.of.Geraldine)
all.baskets

## Repeating vectors
rep(c(0, 0, 7), times = 3)
rep(c(2, 4, 2), each = 3)
rep(c(0, 7), times = c(4,2))
rep(1:3,length.out=7)

# Getting Values in and out of Vectors

## Understanding indexing in R

numbers <- 30:1
numbers

## Extracting values from a vector

numbers[5]
numbers[c(5,11,3)]

indices <- c(5,11,3)
numbers[indices]
numbers[-3]
numbers[-(1:20)]
# numbers[-1:20] # NOT RUN, gives error


## Changing values in a vector


baskets.of.Granny[3] <- 5
baskets.of.Granny

baskets.of.Geraldine[c(2,4)] <- 4
baskets.of.Geraldine

Granny.copy <- baskets.of.Granny

baskets.of.Granny[4] <- 11
baskets.of.Granny

baskets.of.Granny <- Granny.copy
baskets.of.Granny

# Working with Logical Vectors

## Comparing values

baskets.of.Granny > 5
which(baskets.of.Granny > 5)

the.best <- baskets.of.Geraldine < baskets.of.Granny
which(the.best)



## Using logical vectors as indices

baskets.of.Granny[the.best]
x <- c(3, 6, 1, NA, 2)
x[x > 2]
x > 2

## Combining logical statements

min.baskets <- baskets.of.Granny == min(baskets.of.Granny)
max.baskets <- baskets.of.Granny == max(baskets.of.Granny)
min.baskets | max.baskets

x[!is.na(x)]

## Summarizing logical vectors


sum(the.best)
any(the.best)
all(the.best)




# Powering Up Your Math with Vector Functions


## Using arithmetic vector operations

### Summarizing a vector
min(baskets.of.Granny)

max(baskets.of.Granny)
sum(baskets.of.Granny,baskets.of.Geraldine)



x <- c(3,6,2,NA,1)
sum(x)
sum(x,na.rm=TRUE)

### Cumulating operations

cumsum(baskets.of.Granny)
cummax(baskets.of.Geraldine)
cummin(x)

### Calculating differences

diff(baskets.of.Granny)
diff(x)

## Recycling arguments

Granny.pointers <- c(10,2,4,0,4,1,4,2,7,2,1,2)
points <- Granny.pointers * c(2,3)
points
sum(points)

sum(Granny.pointers * c(2,3))



round(diff(baskets.of.Granny) / baskets.of.Granny * 100 )
round(diff(baskets.of.Granny) / baskets.of.Granny[1:5] * 100)

