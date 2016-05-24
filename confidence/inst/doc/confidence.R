## ----ini, echo=FALSE, results='hide', message=FALSE----------------------
library(confidence)
library(knitr)
library(xtable)

## ----echo=FALSE----------------------------------------------------------
opts_chunk$set(
    comment = NA,
    quiet = TRUE,
    progress = FALSE,
    tidy = FALSE,
    cache = FALSE,
    message = FALSE,
    error = TRUE,
    warning = TRUE
)

## ----echo=FALSE----------------------------------------------------------
rversion <- sub(
    pattern = "R *\\(>= *([^)]*)\\).*", 
    replacement = "\\1", 
    x = packageDescription("confidence", fields = "Depends")
)

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  library(confidence)

## ----results='asis', echo=FALSE------------------------------------------
data(metal)
metal$sampdev <- NULL
metal$char <- NULL
metal$comp <- NULL
print(
    xtable(x = metal), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----eval=FALSE, prompt=TRUE---------------------------------------------
#  conf()

## ----eval=FALSE----------------------------------------------------------
#  "my_directory/my_input_file.csv"

## ----eval=FALSE----------------------------------------------------------
#  "my_directory/outputYYYYmmddHHMMSS/output.csv"
#  "my_directory/outputYYYYmmddHHMMSS/output.html"

## ----eval=FALSE, prompt=TRUE---------------------------------------------
#  conf("my_directory/my_input_file.csv")

## ----eval=FALSE, prompt=TRUE---------------------------------------------
#  for (filename in filenames) {
#      conf(filename)
#  }

## ----eval=FALSE, prompt=TRUE---------------------------------------------
#  conf(my_data.frame)

## ----results='asis', echo=FALSE------------------------------------------
data(DCA)
print(
    xtable(x = DCA), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----results='asis', echo=FALSE------------------------------------------
x <- mya(conf_input(x=DCA))
print(
    xtable(x = as.data.frame(x)), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----dev='pdf',fig.width=7, fig.height=4, out.width='0.8\\textwidth', echo=FALSE----
plot(x)

## ----prompt=TRUE---------------------------------------------------------

# load confidence package
# Note: this has to be done only once, at the start of an R session
library(confidence)

# load ecological quality ratio's
data(EQR)

# print these data to the screen
EQR

