## ---- echo=FALSE, results='hide'-----------------------------------------
library("rio")

export(iris, "iris.csv")
export(iris, "iris.rds")
export(iris, "iris.dta")
export(iris, "iris_noext", format = "csv")

## ------------------------------------------------------------------------
library("rio")

x <- import("iris.csv")
y <- import("iris.rds")
z <- import("iris.dta")

# confirm identical
all.equal(x, y, check.attributes = FALSE)
all.equal(x, z, check.attributes = FALSE)

## ------------------------------------------------------------------------
head(import("iris_noext", format = "csv"))

## ---- echo=FALSE, results='hide'-----------------------------------------
unlink("iris.csv")
unlink("iris.rds")
unlink("iris.dta")
unlink("iris_noext")

## ------------------------------------------------------------------------
library("rio")

export(iris, "iris.csv")
export(iris, "iris.rds")
export(iris, "iris.dta")

## ------------------------------------------------------------------------
library("magrittr")
mtcars %>% subset(hp > 100) %>%  aggregate(. ~ cyl + am, data = ., FUN = mean) %>% export(file = "mtcars2.dta")

## ---- echo=FALSE, results='hide'-----------------------------------------
unlink("iris.csv")
unlink("iris.rds")
unlink("iris.dta")
unlink("mtcars2.dta")

## ------------------------------------------------------------------------
# create file to convert
export(iris, "iris.dta")

# convert Stata to SPSS
convert("iris.dta", "iris.sav")

## ------------------------------------------------------------------------
# create an ambiguous file
fwf <- tempfile(fileext = ".fwf")
cat(file = fwf, "123456", "987654", sep = "\n")

# see two ways to read in the file
identical(import(fwf, widths = c(1,2,3)), import(fwf, widths = c(1,-2,3)))

# convert to CSV
convert(fwf, "fwf.csv", in_opts = list(widths = c(1,2,3)))
import("fwf.csv") # check conversion

## ---- echo=FALSE, results='hide'-----------------------------------------
unlink("iris.dta")
unlink("iris.sav")
unlink("fwf.csv")
unlink(fwf)

