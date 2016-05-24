## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(tidy = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  install.packages('formatR', repos = 'http://cran.rstudio.com')
#  #' to install the development version, run
#  #' install.packages('formatR', repos = 'http://yihui.name/xran')

## ------------------------------------------------------------------------
library(formatR)
sessionInfo()

## ----example, eval=FALSE, tidy=FALSE-------------------------------------
#  ## comments are retained;
#  # a comment block will be reflowed if it contains long comments;
#  #' roxygen comments will not be wrapped in any case
#  1+1
#  
#  if(TRUE){
#  x=1  # inline comments
#  }else{
#  x=2;print('Oh no... ask the right bracket to go away!')}
#  1*3 # one space before this comment will become two!
#  2+2+2    # only 'single quotes' are allowed in comments
#  
#  lm(y~x1+x2, data=data.frame(y=rnorm(100),x1=rnorm(100),x2=rnorm(100)))  ### a linear model
#  1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line
#  ## here is a long long long long long long long long long long long long long comment which will be wrapped

## ----example, eval=FALSE, tidy.opts=list(width.cutoff=70)----------------
#  ## comments are retained;
#  # a comment block will be reflowed if it contains long comments;
#  #' roxygen comments will not be wrapped in any case
#  1+1
#  
#  if(TRUE){
#  x=1  # inline comments
#  }else{
#  x=2;print('Oh no... ask the right bracket to go away!')}
#  1*3 # one space before this comment will become two!
#  2+2+2    # only 'single quotes' are allowed in comments
#  
#  lm(y~x1+x2, data=data.frame(y=rnorm(100),x1=rnorm(100),x2=rnorm(100)))  ### a linear model
#  1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line
#  ## here is a long long long long long long long long long long long long long comment which will be wrapped

## ----collapse=TRUE-------------------------------------------------------
library(formatR)
usage(glm, width=70)  # can set arbitrary width here
args(glm)

## ----comment=NA----------------------------------------------------------
set.seed(123)
tidy_eval(text = c("a<-1+1;a  # print the value", "matrix(rnorm(10),5)"))

## ----eval=FALSE----------------------------------------------------------
#  library(formatR)
#  tidy_eval()  # without specifying any arguments, it reads code from clipboard

## ----example, eval=FALSE, echo=6, tidy.opts=list(arrow=TRUE)-------------
#  ## comments are retained;
#  # a comment block will be reflowed if it contains long comments;
#  #' roxygen comments will not be wrapped in any case
#  1+1
#  
#  if(TRUE){
#  x=1  # inline comments
#  }else{
#  x=2;print('Oh no... ask the right bracket to go away!')}
#  1*3 # one space before this comment will become two!
#  2+2+2    # only 'single quotes' are allowed in comments
#  
#  lm(y~x1+x2, data=data.frame(y=rnorm(100),x1=rnorm(100),x2=rnorm(100)))  ### a linear model
#  1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line
#  ## here is a long long long long long long long long long long long long long comment which will be wrapped

## ----example, eval=FALSE, echo=1:6, tidy.opts=list(blank = FALSE)--------
#  ## comments are retained;
#  # a comment block will be reflowed if it contains long comments;
#  #' roxygen comments will not be wrapped in any case
#  1+1
#  
#  if(TRUE){
#  x=1  # inline comments
#  }else{
#  x=2;print('Oh no... ask the right bracket to go away!')}
#  1*3 # one space before this comment will become two!
#  2+2+2    # only 'single quotes' are allowed in comments
#  
#  lm(y~x1+x2, data=data.frame(y=rnorm(100),x1=rnorm(100),x2=rnorm(100)))  ### a linear model
#  1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line
#  ## here is a long long long long long long long long long long long long long comment which will be wrapped

## ----example, eval=FALSE, echo=6, tidy.opts=list(indent = 2)-------------
#  ## comments are retained;
#  # a comment block will be reflowed if it contains long comments;
#  #' roxygen comments will not be wrapped in any case
#  1+1
#  
#  if(TRUE){
#  x=1  # inline comments
#  }else{
#  x=2;print('Oh no... ask the right bracket to go away!')}
#  1*3 # one space before this comment will become two!
#  2+2+2    # only 'single quotes' are allowed in comments
#  
#  lm(y~x1+x2, data=data.frame(y=rnorm(100),x1=rnorm(100),x2=rnorm(100)))  ### a linear model
#  1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line
#  ## here is a long long long long long long long long long long long long long comment which will be wrapped

## ----example, eval=FALSE, echo=6, tidy.opts=list(brace.newline = TRUE)----
#  ## comments are retained;
#  # a comment block will be reflowed if it contains long comments;
#  #' roxygen comments will not be wrapped in any case
#  1+1
#  
#  if(TRUE){
#  x=1  # inline comments
#  }else{
#  x=2;print('Oh no... ask the right bracket to go away!')}
#  1*3 # one space before this comment will become two!
#  2+2+2    # only 'single quotes' are allowed in comments
#  
#  lm(y~x1+x2, data=data.frame(y=rnorm(100),x1=rnorm(100),x2=rnorm(100)))  ### a linear model
#  1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line
#  ## here is a long long long long long long long long long long long long long comment which will be wrapped

## ----example, eval=FALSE, tidy.opts=list(comment = FALSE)----------------
#  ## comments are retained;
#  # a comment block will be reflowed if it contains long comments;
#  #' roxygen comments will not be wrapped in any case
#  1+1
#  
#  if(TRUE){
#  x=1  # inline comments
#  }else{
#  x=2;print('Oh no... ask the right bracket to go away!')}
#  1*3 # one space before this comment will become two!
#  2+2+2    # only 'single quotes' are allowed in comments
#  
#  lm(y~x1+x2, data=data.frame(y=rnorm(100),x1=rnorm(100),x2=rnorm(100)))  ### a linear model
#  1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line
#  ## here is a long long long long long long long long long long long long long comment which will be wrapped

## ----comment-brace, tidy=FALSE, eval=FALSE-------------------------------
#  if (TRUE) {## comments
#  }

## ----comment-brace, eval=FALSE-------------------------------------------
#  if (TRUE) {## comments
#  }

