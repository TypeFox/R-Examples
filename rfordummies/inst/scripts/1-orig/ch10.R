# Chapter 10
# Debugging Your Code

# NOTE : Much code is commented out, as they generate
# errors on purpose. Uncomment the code and run the
# line to see the error and try the debugging out

# Knowing What to Look For

# Reading Errors and Warnings

## Reading error messages

# "a" + 1
# Error in "a" + 1 : non-numeric argument to binary operator

# data.frame(1:10,10:1,)
# Error in data.frame(1:10, 10:1, ) : argument is missing, with no default

## Caring about warnings (or not)

x <- 1:10
y <- if (x < 5 ) 0 else 1

x <- 4
sqrt(x - 5)

plot(1:10, 10:1, color='green')


# Going Bug Hunting

## Calculating the logit

# checks input and does logit calculation
logit <- function(x){
  x <- ifelse(x < 0 | x > 1, "NA", x)
  log(x / (1 - x) )
}
# transforms percentage to number and calls logit
logitpercent <- function(x){
  x <- gsub("%", "", x)
  logit(as.numeric(x))
}



## Knowing where an error comes from

# logitpercent('50%')
# Error in 1 - x : non-numeric argument to binary operator

# traceback()

## Looking inside a function


### Telling R which function to debug

# debug(logit)
# logitpercent('50%')

### Stepping through the function

### Start browsing from within the function

logit <- function(x){
  x <- ifelse(x < 0 | x > 1, "NA", x)
  browser()
  log(x / (1 - x) )
}

# logit(50)

# Generating Your Own Messages

## Creating errors

logit <- function(x){
  if( any(x < 0 | x > 1) ) stop('x not between 0 and 1')
  log(x / (1 - x) )
}


# logitpercent(c('50%','150%'))
# Error in logit(as.numeric(x)/100) : x not between 0 and 1

## Creating warnings
# Function wrapped around for illustrative purposes
# In book only body is given
logit <- function(x){
  x <- ifelse(x < 0 | x > 1, NA, x )
  if( any(is.na(x)) ) warning('x not between 0 and 1')
  log(x / (1 - x) )
}

logitpercent(c('50%','150%'))

# Recognizing the Mistakes You're Sure to Make

## Starting with the wrong data

## Having your data in the wrong format



### Dropping dimensions when you don’t expect it

rowsum.df <- function(x){
  id <- sapply(x,is.numeric)
  rowSums(x[, id])
}

# rowsum.df(sleep)

### Messing up with lists
strsplit('this is a sentence',' ')[2]

strsplit('this is a sentence',' ')

strsplit('this is a sentence',' ')[[1]][2]

customer <- c('Johan Delong','Marie Petit')
namesplit <- strsplit(customer,' ')

paste(namesplit[2],collapse='.')

paste(namesplit[[2]],collapse='.')



### Mixing up factors and numeric vectors

cyl.factor <- as.factor(mtcars$cyl)

median(as.numeric(cyl.factor))

as.numeric(levels(cyl.factor))[cyl.factor]
