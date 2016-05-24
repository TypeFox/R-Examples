## ----ini, echo=FALSE, results='hide', message=FALSE----------------------
library(BEQI2)
library(knitr)
library(xtable)

## ----echo=FALSE----------------------------------------------------------
opts_chunk$set(
    echo = FALSE,
    comment = NA,
    quiet = TRUE,
    progress = FALSE,
    tidy = FALSE,
    cache = FALSE,
    message = FALSE,
    error = TRUE,
    warning = TRUE
)

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  BEQI2dir()

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  BEQI2dir(path = "c:/myprojects/BEQI2/BEQI2_FILES")

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  BEQI2()

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  BEQI2(filename = "c:/myprojects/BEQI2/BEQI2_FILES/settings.json")

## ----echo=FALSE, results='asis'------------------------------------------
cat(
    paste(
        readLines(system.file("extdata", "settings.json", package = "BEQI2")), 
        collapse = "\n"
    )
)

## ----echo=FALSE, results='asis'------------------------------------------
d <- scan(
    file = "./tables/tabbeqi2input.csv", 
    what = character(), 
    sep = ",",
    quiet = TRUE
)
h <- d[1:3]
d <- as.data.frame(matrix(data = d[-(1:3)], ncol = 3, byrow = TRUE))
colnames(d) <- h
print(
    xtable(x = d, align = "llp{50mm}p{50mm}"), 
    include.rownames = FALSE, 
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----echo=FALSE, results='asis'------------------------------------------
d <- readAMBI()
print(
    xtable(x = d[sample.int(n = nrow(d), size = 25), ], align = "llr"), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----echo=FALSE, results='asis'------------------------------------------
filename <- system.file(
    "extdata", "REF-FILES", "BEQI2-Ecotopes.csv", 
    package = "BEQI2"
)
d <- read.csv(file = filename)
print(
    xtable(x = d), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}"),
    sanitize.text.function = function(x) {
        gsub(pattern="_", replacement = "\\\\_", x=x)
    },
    rotate.colnames = TRUE
)

