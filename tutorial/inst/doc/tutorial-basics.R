## ----eval = FALSE--------------------------------------------------------
## library(tutorial)
## tutorial::render("my_doc.Rmd")

## ----eval = FALSE, include = FALSE, ex="play_around", type="sample-code"----
## a <- 2
## b <- 3

## ----eval = FALSE, include = FALSE, ex="create_a", type="pre-exercise-code"----
## # This code is available in the workspace when the session initializes
## b <- 5

## ----eval = FALSE, include = FALSE, ex="create_a", type="sample-code"----
## # Create a variable a, equal to 5
## 
## 
## # Print out a
## 

## ----eval = FALSE, include = FALSE, ex="create_a", type="solution"-------
## # Create a variable a, equal to 5
## a <- 5
## 
## # Print out a
## a

## ----eval = FALSE, include = FALSE, ex="create_a", type="sct"------------
## test_object("a")
## test_output_contains("a", incorrect_msg = "Make sure to print `a`")
## success_msg("Great!")

## ----eval = FALSE, include = FALSE, ex="create_a", type="hint"-----------
## Here is a hint: use `<-` for assignment

