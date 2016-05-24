# Chapter 11 - Getting Help


# Finding Information in the R Help Files

## When you know exactly what you’re looking for

?date

## When you don’t know exactly what you’re looking for

??date


# Searching the Web for Help with R

RSiteSearch("cluster analysis")

install.packages("sos")
library("sos")
findFn("cluster")

# Getting Involved in the R Community

## Using the R mailing lists

## Discussing R on Stack Overflow and Stack Exchange

## Tweeting about R

# Making a Minimal Reproducible Example

dput(cars[1:4, ])

## Creating sample data with random values

set.seed(1)
x <- rnorm(5)
x

cards <- c(1:9, "J", "Q", "K", "A")
suits <- c("Spades", "Diamonds", "Hearts", "Clubs")
deck <- paste(rep(suits, each=13), cards)
set.seed(123)
sample(deck, 7)

set.seed(5)
sample(LETTERS[1:3], 12, replace=TRUE)

set.seed(42)
dat <- data.frame(
   x = sample(1:5),
   y = sample(c("yes", "no"), 5, replace = TRUE)
)
dat

dput(cars[1:4, ])

## Producing minimal code

## Providing the necessary information

sessionInfo()
