## ----setup, echo=FALSE--------------------------------------------------------
library('knitr')
opts_chunk$set(comment=NA, message=FALSE)
options(formatR.blank=FALSE, width=80)

## ----roxygen-ex---------------------------------------------------------------
## the source code of the function `parse_and_save`
ex.file = system.file('examples','parse_and_save.R',package='Rd2roxygen')
cat(readLines(ex.file), sep='\n')

## ----rd-ex--------------------------------------------------------------------
rd.file = system.file('examples','parse_and_save.Rd',package='Rd2roxygen')
cat(readLines(rd.file), sep='\n')

## ----Rd2roxygen---------------------------------------------------------------
library(Rd2roxygen)
formatR::usage(Rd2roxygen)
## e.g. Rd2roxygen('somewhere/to/source/pkg')
## there must be 'man' and 'R' directories under this path

## ----parse-file---------------------------------------------------------------
## we can specify the roxygen comments prefix (#' by default)
options(roxygen.comment = "##' ")
str(info <- parse_file(rd.file))
cat(create_roxygen(info), sep='\n') # parse_and_save() combines these two steps

## ----rab, eval=FALSE, code=formatR::usage(rab)--------------------------------
#  rab(pkg, build = TRUE, build.opts = "--no-manual", install = FALSE, install.opts = if (build) "" else "--with-keep.source", 
#      check = FALSE, check.opts = "--as-cran --no-manual", remove.check = TRUE, reformat = TRUE, 
#      before = NULL, ...)

## ----reformat-ex1, eval=FALSE, tidy=FALSE-------------------------------------
#  ## original code
#  rab=function(pkg,build=TRUE,install=FALSE,
#  check=FALSE,check.opts='',remove.check=TRUE,reformat=TRUE,...){}

## ----reformat-ex2, eval=FALSE, tidy=TRUE--------------------------------------
#  ## the reformatted code; note the spaces and indent
#  rab=function(pkg,build=TRUE,install=FALSE,
#  check=FALSE,check.opts='',remove.check=TRUE,reformat=TRUE,...){}

## ----lyx-doc, eval=FALSE------------------------------------------------------
#  system.file('doc', 'Rd2roxygen.Rmd', package='Rd2roxygen')

## ----write-bib, include=FALSE-------------------------------------------------
# automatically write the bibliography database
write_bib(c('Rd2roxygen', 'roxygen2', 'formatR', 'knitr'), file = 'Rd2roxygen.bib')

