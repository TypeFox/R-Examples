# Chapter 7
# Working in More Dimensions

# Adding a Second Dimension

## Discovering a new dimension

### Creating your first matrix

first.matrix <- matrix(1:12, ncol=4)
first.matrix
matrix(1:12, ncol=4, byrow=TRUE)

### Looking at the properties

str(first.matrix)
dim(first.matrix)
length(first.matrix)
my.array <- array(1:24, dim=c(3,4,2))
baskets.of.Granny <- c(12,4,5,6,9,3)
baskets.of.Geraldine <- c(5,4,2,4,12,9)
baskets.team <- rbind(baskets.of.Granny, baskets.of.Geraldine)

attributes(my.array)
attr(baskets.team,'season') <- '2010-2011'
attr(baskets.team,'season')
attr(baskets.team,'season') <- NULL

## Combining vectors into a matrix

baskets.of.Granny <- c(12,4,5,6,9,3)
baskets.of.Geraldine <- c(5,4,2,4,12,9)
baskets.team <- rbind(baskets.of.Granny, baskets.of.Geraldine)

baskets.team

cbind(1:3, 4:6, matrix(7:12, ncol=2))

# Using the Indices

## Extracting values from a matrix

### Using numeric indices

first.matrix[1:2, 2:3]
first.matrix[2:3,]

### Dropping values using negative indices

first.matrix[-2,-3]

nr <- nrow(first.matrix)
id <- nr*2+2
first.matrix[-id]

first.matrix[-(2 * nrow(first.matrix) + 2)]


### Juggling dimensions

first.matrix[-c(1, 3), ]
first.matrix[2, , drop=FALSE]

## Replacing values in a matrix

first.matrix[3, 2] <- 4
first.matrix

first.matrix[2, ] <- c(1,3)
first.matrix

first.matrix[1:2, 3:4] <- c(8,4,2,1)
first.matrix

# Naming Matrix Rows and Columns

## Changing the row and column names

rownames(baskets.team) <- c('Granny','Geraldine')
rownames(baskets.team)
colnames(baskets.team) <- c('1st','2nd','3th','4th','5th','6th')
baskets.team

colnames(baskets.team)[3] <- '3rd'

baskets.copy <- baskets.team
colnames(baskets.copy) <- NULL
baskets.copy

## Using names as indices

baskets.team[, c("2nd","5th")]

baskets.team["Granny",]

# Calculating with Matrices

## Using standard operations with matrices
first.matrix + 4

second.matrix <- matrix(1:3, nrow=3, ncol=4)

first.matrix + second.matrix

# first.matrix + second.matrix[,1:3] # gives error for illustration
# Error in first.matrix + second.matrix[, 1:3] : non-conformable arrays

first.matrix + 1:3

## Calculating row and column summaries

rowSums(baskets.team)

## Doing matrix arithmetic

### Transposing a matrix

t(first.matrix)

t(1:10)

t(first.matrix[2,])

### Inverting a matrix

square.matrix <- matrix(c(1,0,3,2,2,4,3,2,1),ncol=3)
solve(square.matrix)

### Multiplying two matrices

first.matrix %*% t(second.matrix)

first.matrix %*% 1:4
1:3 %*% first.matrix

# Adding More Dimensions

## Creating an array

### Using the creator functions

my.array <- array(1:24, dim=c(3,4,2))
my.array

### Changing the dimensions of a vector


my.vector <- 1:24
dim(my.vector) <- c(3,4,2)
identical(my.array, my.vector)

## Using dimensions to extract values

my.array[2,3,1]

my.array[, 3, 2, drop=FALSE]


my.array[2, , ]


# Combining Different Types of Values in a Data Frame

## Creating a data frame from a matrix

### Using the function as.data.frame

baskets.df <- as.data.frame(t(baskets.team))

### Looking at the structure of a data frame

baskets.df
str(baskets.df)

### Counting values and variables

nrow(baskets.df)
length(baskets.df)

## Creating a data frame from scratch

### Making a data frame from vectors

employee <- c('John Doe','Peter Gynn','Jolie Hope')
salary <- c(21000, 23400, 26800)
startdate <- as.Date(c('2010-11-1','2008-3-25','2007-3-14'))

employ.data <- data.frame(employee, salary, startdate)

str(employ.data)

### Keeping characters as characters

employ.data <- data.frame(employee, salary, startdate, stringsAsFactors=FALSE)
str(employ.data)

## Naming variables and observations

### Working with variable names

colnames(employ.data)
names(employ.data)

names(employ.data)[3] <- 'firstday'
names(employ.data)

### Naming observations

rownames(employ.data)
rownames(employ.data) <- c('Chef','BigChef','BiggerChef')
employ.data

# Manipulating Values in a Data Frame

## Extracting variables, observations, and values

### Pretending it's a matrix

baskets.df['3rd', 'Geraldine']
baskets.df[, 1]

str(baskets.df[, 1, drop=FALSE])

### Putting your dollar where your data is

baskets.df$Granny

## Adding observations to a data frame

### Adding a single observation

result <- rbind(baskets.df, c(7,4))
result

baskets.df <- rbind(baskets.df,'7th' = c(7,4))
baskets.df

### Adding a series of new observations using rbind

new.baskets <- data.frame(Granny=c(3,8),Geraldine=c(9,4))
rownames(new.baskets) <- c('8th','9th')
baskets.df <- rbind(baskets.df, new.baskets)

### Adding a series of values using indices

baskets.df[c('8th','9th'), ] <- matrix(c(3,8,9,4), ncol=2)
baskets.df[c('8th','9th'), ] <- c(3,8,9,4)

## Adding variables to a data frame

### Adding a single variable

baskets.of.Gabrielle <- c(11,5,6,7,3,12,4,5,9)
baskets.df$Gabrielle <- baskets.of.Gabrielle

head(baskets.df, 4)

### Adding multiple variables using cbind

new.df <- data.frame(
   Gertrude  =  c(3,5,2,1,NA,3,1,1,4),
   Guinevere =  c(6,9,7,3,3,6,2,10,6)
)

head(cbind(baskets.df, new.df), 4)

# Combining Different Objects in a List

## Creating a list

### Creating an unnamed list

baskets.list <- list(baskets.team, '2010-2011')
baskets.list

### Creating a named list

baskets.nlist <- list(scores=baskets.team, season='2010-2011')
baskets.nlist


### Playing with the names of elements

names(baskets.nlist)

### Getting the number of elements

length(baskets.list)

## Extracting elements from lists

### Using [[]]

baskets.list[[1]]
baskets.nlist[['scores']]

### Using []

baskets.list[-1]
baskets.nlist[names(baskets.nlist)=='season']

## Changing the elements in lists

### Changing the value of elements

baskets.nlist[[1]] <- baskets.df
baskets.nlist[['scores']] <- baskets.df
baskets.nlist$scores <- baskets.df

baskets.nlist[1] <- list(baskets.df)

baskets.list[1:2] <- list(baskets.df, '2009-2010')

### Removing elements

baskets.nlist[[1]] <- NULL
baskets.nlist$scores <- NULL
baskets.nlist['scores'] <- NULL

baskets.nlist <- list(scores=baskets.df, season='2010-2011')
baskets.nlist['scores'] <- list(NULL)
baskets.nlist

### Adding extra elements using indices

baskets.nlist$players <- c('Granny','Geraldine')
baskets.nlist[['players']] <- c('Granny','Geraldine')
baskets.nlist['players'] <- list(c('Granny','Geraldine'))

baskets.list[[3]] <- c('Granny','Geraldine')
baskets.list[3] <- list(c('Granny','Geraldine'))

### Combining lists



baskets.list <- list(baskets.team,'2010-2011')
players <- list(rownames(baskets.team))

c(baskets.list, players)

## Reading the output of str() for lists

str(baskets.list)

## Seeing the forest through the trees
