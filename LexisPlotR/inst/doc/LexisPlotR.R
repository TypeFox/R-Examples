## ---- message=FALSE------------------------------------------------------
library(LexisPlotR)

## ------------------------------------------------------------------------
# Plot a Lexis grid from year 1900 to year 1905, representing the ages from 0 to 5
lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)

## ------------------------------------------------------------------------
lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 7)

## ------------------------------------------------------------------------
mylexis <- lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)
mylexis

## ------------------------------------------------------------------------
# Highlight all points that belong to the age of 2
lexis.age(lg = mylexis, age = 2)

## ------------------------------------------------------------------------
# Change the fill colour to "red" and make the layer nearly non-transparent
lexis.age(lg = mylexis, age = 2, fill = "red", alpha = 0.9)

## ------------------------------------------------------------------------
# Highlight the year 1902
lexis.year(lg = mylexis, year = 1902)
# Highlight the cohort 1898
lexis.cohort(lg = mylexis, cohort = 1898)

## ------------------------------------------------------------------------
# Highlight the year 1902, change fill colour to "orange" and increase transparency
lexis.year(lg = mylexis, year = 1902, fill = "orange", alpha = 0.2)
# Highlight the cohort 1898, change fill colour to "grey" and decrease transparency
lexis.cohort(lg = mylexis, cohort = 1898, fill = "grey", alpha = 0.8)

## ------------------------------------------------------------------------
mylexis <- lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)
mylexis <- lexis.age(lg = mylexis, age = 2)
mylexis <- lexis.year(lg = mylexis, year = 1903)
mylexis <- lexis.cohort(lg = mylexis, cohort = 1898)
mylexis

## ------------------------------------------------------------------------
# Define a Lexis grid
mylexis <- lexis.grid(year.start = 1990, year.end = 1995, age.start = 0, age.end = 5)
# Add a life line for an individual born on 1991-09-23
lexis.lifeline(lg = mylexis, entry = "1991-09-23")

## ------------------------------------------------------------------------
lexis.lifeline(lg = mylexis, entry = "1991-09-23", exit = "1994-06-11")

## ------------------------------------------------------------------------
data("lifelines_sample")
str(lifelines_sample)
head(lifelines_sample, 10)

## ------------------------------------------------------------------------
mylexis <- lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)
lexis.lifeline(lg = mylexis, entry = lifelines_sample$entry, exit = lifelines_sample$exit)

## ------------------------------------------------------------------------
lexis.lifeline(lg = mylexis, entry = lifelines_sample$entry, exit = lifelines_sample$exit, lineends = TRUE, colour = "blue", lwd = 1.5, alpha = 0.3)

## ------------------------------------------------------------------------
# Find the path to the sample data
path <- system.file("extdata", "Deaths_lexis_sample.txt", package = "LexisPlotR")
# read the raw data with prepare.hmd()
mydata <- prepare.hmd(path)
# Inspect your data
str(mydata)
summary(mydata[,c("Year", "Age", "Cohort")])

## ---- fig.width=10-------------------------------------------------------
mylexis <- lexis.grid(year.start = 1980, year.end = 1985, age.start = 0, age.end = 5)
# Plot total death counts
lexis.hmd(lg = mylexis, hmd.data = mydata, column = "Total")

## ---- fig.width=10-------------------------------------------------------
mydata$ratioMales <- mydata$Male / mydata$Total
lexis.hmd(lg = mylexis, hmd.data = mydata, column = "ratioMales")

## ------------------------------------------------------------------------
mylexis <- lexis.grid(year.start = 1900, year.end = 1905, age.start = 0, age.end = 5)
# Add a title
mylexis <- mylexis + labs(title = "LexisPlotR")
mylexis
# Change axis labels
mylexis <- mylexis + theme(axis.title = element_text(face = "bold", colour = "red"))
mylexis

