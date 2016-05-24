## ----eval = FALSE, echo = TRUE-------------------------------------------
#  library(GetoptLong)
#  
#  cutoff = 0.05
#  GetoptLong(
#      "number=i", "Number of items, integer, mandatory option",
#      "cutoff=f", "cutoff to filter results, optional, default (0.05)",
#      "verbose",  "print messages"
#  )

## ----eval = FALSE, echo = TRUE, results = 'makeup', highlight = FALSE----
#  ~\> Rscript test.R --number 4 --cutoff 0.01 --verbose
#  ~\> Rscript test.R -n 4 -c 0.01 -v
#  ~\> Rscript test.R -n 4 --verbose

## ----eval = FALSE, echo = TRUE, results = 'makeup', highlight = FALSE----
#  length|size|l=i@

## ----eval = FALSE, echo = TRUE, results = 'makeup', highlight = FALSE----
#  ~\> Rscript foo.R --length 1
#  ~\> Rscript foo.R -len 1
#  ~\> Rscript foo.R --size 1
#  ~\> Rscript foo.R -l 1

## ----eval = FALSE, echo = TRUE-------------------------------------------
#  GetoptLong.options('startingMsg' =
#  'An example to show how to use the packages
#  ')
#  
#  GetoptLong.options('endingMsg' =
#  'Please contact author@gmail.com for comments
#  ')
#  VERSION = "0.0.1"
#  GetoptLong(...)

## ----eval = FALSE, echo = TRUE, results = 'makeup', highlight = FALSE----
#  ~\> Rscript command.R --help
#  An example to show how to use the packages
#  Usage: Rscript test.R [options]
#  
#    --tag integer
#      this is a description of tag which is long long and very long and extremly
#      long...
#  
#    --help
#      Print help message and exit
#  
#    --version
#      Print version information and exit
#  
#  Please contact author@gmail.com for comments

## ----eval = FALSE, echo = TRUE, results = 'makeup', highlight = FALSE----
#  ~\> Rscript command.R --version
#  0.0.1

## ----eval = FALSE, echo = TRUE, results = 'hide'-------------------------
#  GetoptLong.options("config" = "bundling")
#  GetoptLong.options("config" = c("no_ignore_case", "bundling"))

## ----eval = FALSE, echo = TRUE, results = 'makeup', highlight = FALSE----
#  -a -b -c  -abc
#  -s 24 -s24 -s=24

## ----eval = FALSE, echo = TRUE, results = 'makeup', highlight = FALSE----
#  ~\> Rscript test.R -a -b -c -- /your/perl/bin/perl

