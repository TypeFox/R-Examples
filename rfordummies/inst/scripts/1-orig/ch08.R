# Chapter 8
# Putting the Fun in Functions

# Moving from Scripts to Functions

## Making the script

x <- c(0.458, 1.6653, 0.83112)
percent <- round(x * 100, digits = 1)
result <- paste(percent, "%", sep = "")
print(result)

# source('pastePercent.R') # Only after saving

## Transforming the script

addPercent <- function(x){
  percent <- round(x * 100, digits = 1)
  result <- paste(percent, "%", sep = "")
  return(result)
}

## Using the function

ls()

### Formatting the numbers

new.numbers <- c(0.8223, 0.02487, 1.62, 0.4)
addPercent(new.numbers)

### Playing with function objects

ppaste <- addPercent
ppaste

## Reducing the number of lines

### Returning values by default

# AddPercent function without last return - not written in book
addPercent <- function(x){
  percent <- round(x * 100, digits = 1)
  result <- paste(percent, "%", sep = "")
}

print( addPercent(new.numbers) )

addPercent <- function(x){
  percent <- round(x * 100, digits = 1)
  paste(percent, "%", sep = "")
}

addPercent <- function(x){
  if( !is.numeric(x) ) return(NULL)
  percent <- round(x * 100, digits = 1)
  paste(percent, "%", sep = "")
}

### Breaking the walls

odds <- function(x) x / (1-x)

odds(0.8)

addPercent <- function(x) paste(round(x * 100, digits = 1), "%", sep = "")

# Using Arguments the Smart Way

## Adding more arguments

percentages <- c(58.23, 120.4, 33)
addPercent(percentages/100)

### Adding the mult argument

addPercent <- function(x, mult){
  percent <- round(x * mult, digits = 1)
  paste(percent, "%", sep = "")
}

addPercent(percentages, mult = 1)

### Adding a default value

# addPercent(new.numbers) # Gives error for illustrative purposes
# Error in x * mult : 'mult' is missing

addPercent <- function(x, mult = 100){
  percent <- round(x * mult, digits = 1)
  paste(percent, "%", sep = "")
}

addPercent(new.numbers)

addPercent(percentages, 1)

## Conjuring tricks with dots

addPercent <- function(x, mult = 100, ...){
  percent <- round(x * mult, ...)
  paste(percent, "%", sep = "")
}

addPercent(new.numbers, digits = 2)
addPercent(new.numbers)


addPercent <- function(x, mult = 100, digits = 1){
  percent <- round(x * mult, digits = digits)
  paste(percent, "%", sep = "")
}

## Using functions as arguments

### Applying different ways of rounding

addPercent <- function(x, mult = 100, FUN = round, ...){
  percent <- FUN(x * mult, ...)
  paste(percent, "%", sep = "")
}

addPercent(new.numbers, FUN = signif, digits = 3)

### Using anonymous functions

profits <- c(2100, 1430, 3580, 5230)
rel.profit <- function(x) round(x / sum(x) * 100)
addPercent(profits,
                FUN = function(x) round(x / sum(x) * 100) )

addPercent(profits / sum(profits))

# Coping with Scoping



## Crossing the borders

### Creating a test case

x <- 1:5
test <- function(x){
  cat("This is x:", x, "\n")
  rm(x)
  cat("This is x after removing it:",x,"\n")
}

test(5:1)

### Searching the path

## Using internal functions

calculate.eff <- function(x, y, control){
  min.base <- function(z) z - mean(control)
  min.base(x) / min.base(y)
}


half <- c(2.23, 3.23, 1.48)
full <- c(4.85, 4.95, 4.12)
nothing <- c(0.14, 0.18, 0.56, 0.23)
calculate.eff(half, full, nothing)

# Dispatching to a Method

## Finding the methods behind the function

print

### Using methods with UseMethod

small.one <- data.frame(a = 1:2, b = 2:1)
print.data.frame(small.one)

### Using default methods

print.default(small.one)

## Doing it yourself

### Adapting the addPercent function

addPercent.character <- function(x){
  paste(x,"%",sep="")
}

# Not written out in the book - needed for rest code # 
addPercent.numeric <- function(x, mult = 100, FUN = round, ...){
  percent <- FUN(x * mult, ...)
  paste(percent, "%", sep = "")
}

addPercent <- function(x,...){
  UseMethod("addPercent")
}

addPercent(new.numbers, FUN = floor)

addPercent(letters[1:6])

# Adding a default function

# addPercent(small.one) # Gives error on purpose
# Error in UseMethod("addPercent") :
#  no applicable method for 'addPercent' applied to an object of class "data.frame"

addPercent.default <- function(x){
  cat('You should try a numeric or character vector.\n')
}