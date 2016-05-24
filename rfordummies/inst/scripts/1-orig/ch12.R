# Chapter 12
# Getting Data into and out of R

# NOTE : Most of the code depends on actions, directories
# and the presence of files. Code that isn't runnable is
# commented out.

# Getting Data into R

## Entering data in the R text editor

elements <- data.frame()
# elements <- edit(elements)

# print(elements)

## Using the Clipboard to copy and paste
# Reminder : This only works on Windows

# x <- readClipboard()
# x
# x <- readClipboard()
# x
# x <- read.table(file = "clipboard", sep = "\t", header=TRUE)
# x

## Reading data in CSV files

### Using read.csv() to import data

# elements <- read.csv(file.path("f:", "elements.csv"))
# str(elements)
# elements <- read.csv(file.path("f:", "elements.csv"), stringsAsFactors=FALSE)
# str(elements)

### Using read.table() to import tabular data in text files

## Reading data from Excel
# install.packages("XLConnect")
# library("XLConnect")
# excel.file <- file.path("~/Elements.xlsx")

# elements <- readWorksheetFromFile(excel.file, sheet=1)
# elements <- readWorksheetFromFile(excel.file, sheet="Sheet1")

## Working with other data types

# library(foreign)
# read.spss(file="location/of/myfile")

# Getting Your Data out of R


# writeClipboard(names(iris))

# write.table(head(iris), file="clipboard", sep="\t", row.names=FALSE)

# Working with Files and Folders

## Understanding the working directory
getwd()

# setwd("F:/git/roxygen2")
# getwd()
# setwd("F:\\git\\stringr")
# getwd()

file.path("f:", "git", "surveyor")

# setwd(file.path("F:", "git", "roxygen2"))
# getwd()

file.path("F:", "git", "roxygen2", "roxygen2", "README.md" )

## Manipulating files

# list.files(file.path("F:", "git", "roxygen2"))
my.file <- tempfile()
my.file
write.csv(iris, file=my.file)
list.files(tempdir())

file.iris <- read.csv(my.file)

file.remove(my.file)
list.files(tempdir())


