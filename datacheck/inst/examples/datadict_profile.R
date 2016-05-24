library(stringr)
# Get example data files
atable <- system.file("examples/db.csv", package = "datacheck")
arule <- system.file("examples/rules1.R", package = "datacheck")
aloctn <- system.file("examples/location.csv", package = "datacheck")  # for use in is.oneOf

ctable <- basename(atable)
crule <- basename(arule)
cloctn <- basename(aloctn)

cwd <- tempdir()
owd <- getwd()
setwd(cwd)

file.copy(atable, ctable)
file.copy(arule, crule)
file.copy(aloctn, cloctn)

at <- read.csv(ctable, stringsAsFactors = FALSE)
ad <- read_rules(crule)

db <- datadict_profile(at, ad)

is_datadict_profile(db) == TRUE

db

setwd(owd)

