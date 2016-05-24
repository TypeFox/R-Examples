
# Get example data file with some errors in it
atbler <- system.file("examples/db-err.csv", package = "datacheck")
arule <- system.file("examples/rules1.R", package = "datacheck")

at <- read.csv(atbler, stringsAsFactors = FALSE)
ad <- read_rules(arule)

db_e <- datadict_profile(at, ad)

has_rule_errors(db_e) == TRUE 
