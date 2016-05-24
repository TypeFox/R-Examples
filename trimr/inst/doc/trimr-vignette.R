## ------------------------------------------------------------------------
# load the trimr package
library(trimr)

# activate the data
data(exampleData)

# look at the top of the data
head(exampleData)

## ------------------------------------------------------------------------
# perform the trimming
trimmedData <- absoluteRT(data = exampleData, minRT = 150, maxRT = 2000, 
                          digits = 0)

# look at the top of the data
head(trimmedData)

## ------------------------------------------------------------------------
# perform the trimming
trimmedData <- absoluteRT(data = exampleData, minRT = 150, maxRT = 2000, 
                          returnType = "raw", digits = 0)

# look at the top of the data
head(trimmedData)

## ------------------------------------------------------------------------
# trim the data
trimmedData <- sdTrim(data = exampleData, minRT = 150, sd = 3, 
                      perCondition = FALSE, perParticipant = TRUE, 
                      returnType = "mean", digits = 0)

# look at the top of the data
head(trimmedData)

## ------------------------------------------------------------------------
# trim the data
trimmedData <- sdTrim(data = exampleData, minRT = 150, sd = 3, 
                      perCondition = TRUE, perParticipant = TRUE, 
                      returnType = "mean", digits = 0)

# look at the top of the data
head(trimmedData)

## ------------------------------------------------------------------------
# load the data
data(linearInterpolation)

# show the first 20 rows (there are 100 in total)
linearInterpolation[1:20, ]

## ------------------------------------------------------------------------
# trim the data
trimmedData <- nonRecursive(data = exampleData, minRT = 150, digits = 0)

# see the top of the data
head(trimmedData)

## ------------------------------------------------------------------------
# trim the data
trimmedData <- modifiedRecursive(data = exampleData, minRT = 150, digits = 0)

# see the top of the data
head(trimmedData)

## ------------------------------------------------------------------------
# trim the data
trimmedData <- hybridRecursive(data = exampleData, minRT = 150, digits = 0)

# see the top of the data
head(trimmedData)

## ------------------------------------------------------------------------
# get the example data that ships with trimr
data(exampleData)

# pass it to a new variable
newData <- exampleData

# add a column called "taskSequence" that logs whether the current task was a 
# repetition or a switch trial (which is currently coded in the "condition")
# column
newData$taskSequence <- newData$condition

# add a column called "reward" that logs whether the participant received a 
# reward or not. Fill it with random entries, just for example. This uses R's
# "sample" function
newData$reward <- sample(c("Reward", "NoReward"), nrow(newData), 
                         replace = TRUE)

# delete the "condition" column
newData <- subset(newData, select = -condition)

# now let's look at our new data
head(newData)

## ------------------------------------------------------------------------
# add a new column called "condition", and fill it with information from both 
# columns that code for our factors
newData$condition <- paste(newData$taskSequence, "_", newData$reward, sep = "")

# let's again look at the data
head(newData)

## ------------------------------------------------------------------------
# trim the data
trimmedData <- sdTrim(newData, minRT = 150, sd = 2.5)

# check it worked
head(trimmedData)


