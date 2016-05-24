library("R.filesets")

message("*** TabularTextFileSet")

# Setup a file set consisting of all *.dat tab-delimited files
# in a particular directory
pathA <- system.file("exData", "dataSetA,original", package="R.filesets")
ds <- TabularTextFileSet$byPath(pathA, pattern="[.]dat$")
print(ds)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract one column with a particular name (one per file)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read column 'y' and a subset of the rows from each of the
# tab-delimited files and combine into a matrix
rows <- c(3:5, 8, 2)
data <- extractMatrix(ds, column="y", colClass="integer", rows=rows, drop=TRUE)
print(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read data frames from each of the files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataList <- lapply(ds, FUN=readDataFrame)
print(dataList)

rows <- c(3:5, 8, 2)
dataList <- lapply(ds, FUN=readDataFrame, rows=rows)
print(dataList)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read common columns and stack into one data frame
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
colNames <- Reduce(intersect, lapply(ds, getColumnNames))
cat("Common column names:\n")
print(colNames)

# Read the *common* columns "as is" (hence 'NA')
colClasses <- rep(NA, times=length(colNames))
names(colClasses) <- colNames
cat("Column class patterns:\n")
print(colClasses)

data <- readDataFrame(ds, colClasses=colClasses, verbose=TRUE)
print(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Translate column names on the fly
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lapply(ds, FUN=setColumnNamesTranslator, function(names, ...) toupper(names))
data <- readDataFrame(ds, colClasses=c("(X|Y)"="integer", "CHAR"="character"))
print(data)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ADVANCED: Translation of fullnames
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fnts <- TabularTextFileSet$byPath(getPath(ds), pattern=",fullnames[.]txt$")
appendFullNamesTranslator(ds, as.list(fnts))

cat("Default fullnames:\n")
print(head(getFullNames(ds, translate=FALSE)))
cat("Translated fullnames:\n")
print(head(getFullNames(ds)))

cat("Default fullnames:\n")
print(getFullNames(ds, translate=FALSE))
cat("Translated fullnames:\n")
print(getFullNames(ds))
