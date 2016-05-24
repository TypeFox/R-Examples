# Chapter 3 - The Fundamentals of R

# Using the Full Power of Functions

## Vectorizing your functions

baskets.of.Granny <- c(12,4,4,6,9,3)
baskets.of.Granny
sum(baskets.of.Granny)

firstnames <- c("Joris", "Carolien", "Koen")
lastname <- "Meys"
paste(firstnames,lastname)

authors <- c("Andrie","Joris")
lastnames <- c("de Vries","Meys")
paste(authors,lastnames)

## Putting the argument in a function

# print() ### This line of code leads to deliberate error for illustration
print(x = "Isn't this fun?")

print(digits=4, x = 11/7)

# Making history

savehistory(file = "Chapter3.Rhistory")
loadhistory("Chapter3.Rhistory")

# Keeping Your Code Readable

## Following naming conventions

## Choosing a clear name

paste <- paste("This gets","confusing")
paste
paste("Don't","you","think?")

## Choosing a naming style

## Structuring your code

baskets.of.Geraldine <- c(5,3,2,2,12,9)
Intro <- "It is amazing! The All Star Grannies scored
a total of"

Outro <- "baskets in the last six games!"

Total.baskets <- baskets.of.Granny +
               baskets.of.Geraldine

Text <- paste(Intro,
              sum(Total.baskets),
              Outro)
cat(Text)
Text

cat('If you doubt whether it works,
+ just try it out.')

## Adding comments

# The All Star Grannies do it again!
baskets.of.Granny <- c(12,4,4,6,9,3) # Granny rules
sum(baskets.of.Granny) # total number of points


# Getting from Base R to More

## Finding packages

## Installing packages

install.packages("fortunes")

library("fortunes")
fortune("This is R")
fortune(161)
detach(package:fortunes)

