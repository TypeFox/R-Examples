##---
##title: "r2rmd function usage example"
##author: "Peter Solymos"
##date: "May 13, 2015"
##output: pdf_document
##---

### Rules

## * leading `##` is treated as *non-code*
## * leading `#` followed by other than `#` is *code comment*
## * leading `#` after whitespace is *code comment*
## * newline is *code* when preceded and followed by code

## The leading double hash `##` trimmed for comment lines.
## R markdown chunk start/end stuff is added for code chunks.
## The argument `extra` adds chunk arguments, e.g. `extra=', eval=FALSE'` etc.
## The function returns vector if `out=NULL`.

### Trivial example

## Addition
a <- 1
b <- a + 1 # comment
## Function
# another comment
f <- function(x) x + a
f (b)
## That's all folks!

