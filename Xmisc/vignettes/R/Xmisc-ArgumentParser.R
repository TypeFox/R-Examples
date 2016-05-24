#!/usr/bin/env Rscript

## ************************************************************************
## This is an executable R script to illustrate `ArgumentParser'
## in CRAN/R package `Xmisc'.
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Mon Aug 11 15:28:28 EDT 2014 -0400 (Week 32)
## 
## 
## Reference:
## CRAN/R package `Xmisc'
## http://CRAN.R-project.org/package=Xmisc
##
## Get help:
## Rscript Xmisc-argumentparser.R -h
## Runme: 
## Rscript Xmisc-argumentparser.R --a_name=${USER} --a_int=2 --a_num=3.6 --if.test=TRUE
## 
## ************************************************************************

require(methods)
require(Xmisc)



PARSEME <- function(){
  parser <- ArgumentParser$new()

  parser$add_usage('Xmisc-argumentparser.R [options]')
  parser$add_description(
    'An executable R script parsing arguments from Unix-like command line.')

  parser$add_argument(
    '--h',type='logical',
    action='store_true',
    help='Print the help page'
  )

  parser$add_argument(
    '--help',type='logical',
    action='store_true',
    help='Print the help page'
  )

  parser$add_argument(
    '--a_name',type='character',
    required=TRUE,
    help='A a_name.'
  )
  
  parser$add_argument(
    '--a_int',type='integer',
    default=1,    
    help='A integer.'
  )

  parser$add_argument(
    '--a_num',type='numeric',
    default=1,    
    help='A number.'
  )

  parser$add_argument(
    '--if.test',type='logical',
    default=FALSE,
    help='Whether it is a test?!'
  )
  
  return(parser)
}


hello <- function(){
  message('Hello, ',a_name,'!')
  message('The integer is ',a_int,'.')
  message('The number is ',a_num,'.')
  message('')
  flag <- if (if.test) "" else "not "
  message('(This is ',flag,'a test).')  
  message('')
}

info <- function(){
  message('class(a_name): ',class(a_name))
  message('class(a_int): ',class(a_int))
  message('class(a_num): ',class(a_num))
  message('class(if.test): ',class(if.test))  
  message('')
}


main <- function(){
  parser <- PARSEME()
  parser$helpme()
  
  hello()
  info()
  
  invisible()
}


main()
