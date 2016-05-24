OneTrueLove <- read.file("OneTrueLove.csv")
head(OneTrueLove)
tally(Response~Gender, format = "count", margins = TRUE, data = OneTrueLove)
prop(Response~Gender, data = OneTrueLove)
diff(prop(Response~Gender, data = OneTrueLove))

