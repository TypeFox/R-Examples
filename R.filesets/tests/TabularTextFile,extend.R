library("R.filesets")

message("*** Extending TabularTextFile")

pathA <- system.file("exData", "dataSetA,original", package="R.filesets")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plain tab-delimited file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tsf <- TabularTextFile("fileB,other,tags.dat", path=pathA)
print(tsf)

# Read all data
data <- readDataFrame(tsf)
str(data)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Custom tab-delimited file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!is.element("covr", loadedNamespaces())) {
  setConstructorS3("MyTabularTextFile", function(...) {
    extend(TabularTextFile(...), "MyTabularTextFile")
  })

  setMethodS3("getDefaultColumnClassPatterns", "MyTabularTextFile", function(...) {
    c("*"="NULL", "(x|y)"="integer")
  })

  tsfx <- MyTabularTextFile("fileB,other,tags.dat", path=pathA)
  print(tsfx)

  rargs <- getReadArguments(tsfx)
  print(rargs$colClasses)
  stopifnot(rargs$colClasses["x"] == "integer",
            rargs$colClasses["y"] == "integer",
            rargs$colClasses["fac"] == "NULL",
            rargs$colClasses["char"] == "NULL")

  datax <- readDataFrame(tsfx)
  str(datax)
}

