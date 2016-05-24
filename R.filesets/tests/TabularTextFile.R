library("R.filesets")

message("*** TabularTextFile")

pathA <- system.file("exData", "dataSetA,original", package="R.filesets")
pathB <- system.file("exData", "dataSetB", package="R.filesets")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# File #1 - regular tab-delimited file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
db <- TabularTextFile("fileB,other,tags.dat", path=pathA)
print(db)

# Read all data
data <- readDataFrame(db, verbose=TRUE)
print(data)

# Read columns
dataC <- readColumns(db, verbose=TRUE)
print(dataC)

# Extract a particular column by its name
dataY <- extractMatrix(db, column="y", colClass="integer")

# Validate
stopifnot(identical(dataY[,1], data$y))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# File #2 - with header comments
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
db <- TabularTextFile("fileA,20100112.dat", path=pathA)
print(db)

# Read all data
data <- readDataFrame(db)
print(data)

# Read columns 'x', 'y', and 'char'
data <- readDataFrame(db, colClasses=c("(x|y)"="integer", "char"="character"))
print(data)

# Translate column names on the fly
db <- setColumnNamesTranslator(db, function(names, ...) toupper(names))
data <- readDataFrame(db, colClasses=c("(X|Y)"="integer", "CHAR"="character"))
print(data)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# File #3 - column names in header comments
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
db <- TabularTextFile("fileE,headerArgs.dat", path=pathA)
print(db)

# Read all data
data <- readDataFrame(db)
print(data)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# File #4 - with neither column names nor header comments
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
db <- TabularTextFile("fileF,noHeader.dat", path=pathB, columnNames=FALSE)
print(db)

# Read all data
data <- readDataFrame(db)
print(data)
str(data)

# Use column classes
colClasses <- rep(NA_character_, times=nbrOfColumns(db))
colClasses[length(colClasses)] <- "NULL"
data <- readDataFrame(db, colClasses=colClasses)
print(data)
str(data)

# Sanity check
stopifnot(ncol(data) == nbrOfColumns(db) - 1L)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# File #5 - with and without newline for the last line
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
df1 <- TabularTextFile("fileG,EOL.txt", path=pathB)
print(df1)
data1 <- readDataFrame(df1)

df2 <- TabularTextFile("fileG,noEOL.txt", path=pathB)
print(df2)
data2 <- readDataFrame(df2)

# Sanity checks
stopifnot(identical(data2, data1))

