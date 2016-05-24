# Chapter 9
# Controlling the Logical Flow

#Making Choices with if Statements

priceCalculator <- function(hours, pph=40){
    net.price <- hours * pph
    round(net.price)
}

priceCalculator <- function(hours, pph=40){
    net.price <- hours * pph
    if(hours > 100) {
      net.price <- net.price * 0.9
    }
    round(net.price)
}

priceCalculator(hours = 55)
priceCalculator(hours = 110)

priceCalculator <- function(hours, pph=40){
    net.price <- hours * pph
    if(hours > 100) net.price <- net.price * 0.9
    round(net.price)
}

?'if'
?"if"
?`if`

## Doing Something Else with an if...else Statement

priceCalculator <- function(hours, pph=40, public=TRUE){
    net.price <- hours * pph
    if(hours > 100) net.price <- net.price * 0.9
    if(public) {
      tot.price <- net.price * 1.06
    } else {
      tot.price <- net.price * 1.12
    }
    round(tot.price)
}

priceCalculator(25,public=TRUE)
priceCalculator(25,public=FALSE)

priceCalculator <- function(hours, pph=40, public=TRUE){
    net.price <- hours * pph
    if(hours > 100) net.price <- net.price * 0.9
    if(public) tot.price <- net.price * 1.06 else
               tot.price <- net.price * 1.12
    round(tot.price)
}

priceCalculator <- function(hours, pph=40, public=TRUE){
    net.price <- hours * pph
    if(hours > 100) net.price <- net.price * 0.9
    tot.price <- net.price * if(public) 1.06 else 1.12
    round(tot.price)
}

# Vectorizing Choices

## Looking at the problem

priceCalculator(c(25,110))
priceCalculator(110)
c(25, 110) > 100

## Choosing based on a logical vector

### Understanding how it works

ifelse(c(1,3) < 2.5 , 1:2 , 3:4)

### Trying it out

my.hours <- c(25,110)
my.hours * 40 * ifelse(my.hours > 100, 0.9, 1)

### Adapting the function

priceCalculator <- function(hours,pph=40,public){
    net.price <- hours * pph
    net.price <- net.price * ifelse(hours > 100 , 0.9, 1)
    tot.price <- net.price * ifelse(public, 1.06, 1.12)
    round(tot.price)
}

clients <- data.frame(
  hours = c(25, 110, 125, 40),
  public = c(TRUE,TRUE,FALSE,FALSE)
)

with(clients, priceCalculator(hours, public = public))

# Making Multiple Choices

## Chaining if...else statements


# Code example # NOT run
#if(client=='private'){
#  tot.price <- net.price * 1.12      # 12% VAT
#} else {
#  if(client=='public'){
#    tot.price <- net.price * 1.06    # 6% VAT
#  } else {
#    tot.price <- net.price * 1    # 0% VAT
#  }
#}

# Code example # NOT run
#if(client=='private'){
#    tot.price <- net.price * 1.12
#} else if(client=='public'){
#    tot.price <- net.price * 1.06
#} else {
#    tot.price <- net.price
#}

# Code example # NOT run
#VAT <- ifelse(client=='private', 1.12,
#          ifelse(client == 'public', 1.06, 1)
#       )
#tot.price <- net.price * VAT
#

## Switching between possibilities

### Making choices with switch

# Code example # NOT run
# VAT <- switch(client, private=1.12, public=1.06, abroad=1)



### Using default values in switch

# Code example # NOT run
# VAT <- switch(client, private=1.12, public=1.06, 1)

client <- 'other'
switch(client, private=1.12, public=1.06, 1)


# Looping Through Values

## Constructing a for loop

## Calculating values in a for loop

### Using the values of the vector

priceCalculator <- function(hours, pph=40, client){
    net.price <- hours * pph *
                   ifelse(hours > 100, 0.9, 1)

    VAT <- numeric(0)
    for(i in client){
      VAT <- c(VAT,switch(i, private=1.12, public=1.06, 1))
    }

    tot.price <- net.price * VAT
    round(tot.price)
}


clients$type <- c('public','abroad','private','abroad')
priceCalculator(clients$hours, client=clients$type)

### Using loops and indices

nclient <- length(client)
VAT <- numeric(nclient)
for(i in seq_along(client)){
  VAT[i] <- switch(client[i], private=1.12, public=1.06, 1)
}
VAT

# Looping without Loops: Meeting the Apply Family

songline <- 'Get out of my dreams...'
for(songline in 1:5) print('...Get into my car!')

songline

## Looking at the family features

## Meeting three of the members


## Applying functions on rows and columns

### Counting birds

counts <- matrix(c(3,2,4,6,5,1,8,6,1), ncol=3)
colnames(counts) <- c('sparrow','dove','crow')
counts

apply(counts, 2, max)

### Adding extra arguments

counts[2, 2] <- NA
apply(counts,2,max)
apply(counts, 2, max, na.rm=TRUE)

## Applying functions to listlike objects

### Applying a function to a vector

#### Using switch on vectors

sapply(c('a','b'), switch, a='Hello', b='Goodbye')

#### Replacing a complete for loop with a single statement

priceCalculator <- function(hours, pph=40, client){
  net.price <- hours * pph * ifelse(hours > 100, 0.9, 1)

  VAT <- sapply(client, switch, private=1.12, public=1.06, 1)

  tot.price <- net.price * VAT
  round(tot.price)
}

### Applying a function to a data frame

sapply(clients,class)

### Simplifying results (or not) with sapply

sapply(clients, unique)

### Getting lists using lapply

sapply(clients[c(1,3), ], unique)

lapply(clients[c(1,3), ], unique)
