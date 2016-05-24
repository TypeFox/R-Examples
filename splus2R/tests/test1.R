library("splus2R")

# Basic test
testFile <- "do.test.t"
do.test(testFile)

# Test verbose=TRUE
do.test(testFile, verbose=TRUE)

# Test strict=TRUE
try(
  do.test(testFile, strict=TRUE)
) # Need try() to prevent R CMD check problems

# Test check argument
do.test(testFile, check = Quote(print(".......")))

# Test use of a connection, as opposed to a file.
testConnection <- file(testFile, open="r")
do.test(testConnection)



#--------------------------------------------------
# Test removal of objects

# Test effect of creating but not removing objects ("x")
testFile2 <- "createObjects.t"
do.test(testFile2)
x
# x should exist in the working environment
exists("x", envir = .GlobalEnv)
if(exists("x", envir = .GlobalEnv))
  rm(x, envir = .GlobalEnv)


do.test(testFile2, local=TRUE)
(!exists("x", envir = .GlobalEnv))
# x should not exist in the working environment



# Test: (1) removal of objects (2) whether objects created are used
x <- 55 # set to 5 in the function; if 5 then tests should pass
testFile3 <- "useObjects.t"
do.test(testFile3)
(!exists("x", envir = .GlobalEnv))
# x should not exist in the working environment

x <- 55 # set to 5 in the function; if 5 then tests should pass
do.test(testFile3, local=TRUE)
all.equal(x, 55)  # is 55 in R
exists("x", envir = .GlobalEnv) && identical(x, 55)
if(exists("x", envir = .GlobalEnv))
  rm(x, envir = .GlobalEnv)


#--------------------------------------------------
# Test auxiliary functions

# Test the use of allTrue
testFile4 <- "allTrue.t"
do.test(testFile4)


# Test the use of all.equal.excluding
testFile5 <- "all.equal.excluding.t"
do.test(testFile5)
