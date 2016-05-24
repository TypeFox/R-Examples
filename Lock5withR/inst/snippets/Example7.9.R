OneTrueLove <- read.file("OneTrueLove.csv")
tally(~Response, format = "proportion", data = OneTrueLove)
tally(~ Response + Gender, format = "proportion", margin = TRUE, data = OneTrueLove)

