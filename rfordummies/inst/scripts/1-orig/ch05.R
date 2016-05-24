# Chapter 5 - Getting Started with Reading and Writing

# Using Character Vectors for Text Data

## Assigning a value to a character vector

x <- "Hello world!"
is.character(x)
length(x)
nchar(x)

## Creating a character vector with more than one element

x <- c("Hello", "world!")
length(x)
nchar(x)

## Extracting a subset of a vector

letters
LETTERS
letters[10]
LETTERS[24:26]
tail(LETTERS, 5)
head(letters, 10)

## Naming the values in your vectors

### Looking at how named vectors work

str(islands)
islands[c("Asia", "Africa", "Antarctica")]
names(islands)[1:9]
names(sort(islands, decreasing=TRUE)[1:6])

## Creating and assigning named vectors

month.days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
names(month.days) <- month.name
month.days
names(month.days[month.days==31])

# Manipulating Text

## String theory: Combining and splitting strings

### Splitting text

pangram <- "The quick brown fox jumps over the lazy dog"
pangram
strsplit(pangram, " ")

words <- strsplit(pangram, " ")[[1]]
words

### Changing text case

unique(tolower(words))
toupper(words[c(4, 9)])
tolower("Some TEXT in Mixed CASE")

### Concatenating text

paste("The", "quick", "brown", "fox")
paste(c("The", "quick", "brown", "fox"))
paste(words, collapse=" ")
paste(words, collapse="_")
paste(LETTERS[1:5], 1:5, sep="_", collapse="---")
paste("Sample", 1:5)
paste(c("A", "B"), c(1, 2, 3, 4), sep="-")
paste(c("A"), c(1, 2, 3, 4, 5), sep="-")

## Sorting text

sort(letters, decreasing=TRUE)
sort(words)

## Finding text inside text

### Searching for individual words

head(state.name)

### Searching by position

head(substr(state.name, start=3, stop=6))

### Searching by pattern

grep("New", state.name)
state.name[29]
state.name[grep("New", state.name)]
state.name[grep("new", state.name)]

### Searching for multiple words

state.name[grep(" ", state.name)]
state.name[grep("East", state.name)]

## Substituting text


gsub("cheap", "sheep's", "A wolf in cheap clothing")
x <- c("file_a.csv", "file_b.csv", "file_c.csv")
y <- gsub("file_", "", x)
y
gsub(".csv", "", y)


#### Extending text functionality with stringr

install.packages("stringr")
library(stringr)


## Revving up with regular expressions

rwords <- c("bach", "back", "beech", "beach", "black")
grep("beach|beech", rwords)
rwords[grep("beach|beech", rwords)]
rwords[grep("be(a|e)ch", rwords)]
rwords[grep("b(e*|a*)ch", rwords)]


# Factoring in Factors

## Creating a factor

directions <- c("North", "East", "South", "South")
factor(directions)
factor(directions, levels= c("North", "East", "South", "West"))
factor(directions, levels= c("North", "East", "South", "West"), labels=c("N", "E", "S", "W"))

## Converting a factor

directions <- c("North", "East", "South", "South")
directions.factor <- factor(directions)
directions.factor
as.character(directions.factor)
as.numeric(directions.factor)

numbers <- factor(c(9, 8, 10, 8, 9))
as.character(numbers)
as.numeric(numbers)
as.numeric(as.character(numbers))

## Looking at levels

str(state.region)
levels(state.region)
levels(state.region) <- c("NE", "S", "NC", "W")
head(state.region)
nlevels(state.region)
length(levels(state.region))
levels(state.region)[2:3]

## Distinguishing data types

head(state.region)
table(state.region)
state.region

## Working with ordered factors

status <- c("Lo", "Hi", "Med", "Med", "Hi")
ordered.status <- factor(status, levels=c("Lo", "Med", "Hi"), ordered=TRUE)
ordered.status
table(status)
table(ordered.status)


