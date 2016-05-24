################################################################################
##
## $Id: df.category.mean.test.R 1625 2010-02-18 19:44:29Z enos $
##
## Tests df.category.mean function
##
################################################################################

library(portfolio)

load("df.category.mean.test.RData")

## save(truth, file = "df.category.mean.test.RData", compress = TRUE)

## Corner cases: data.frame with 0 rows and a data.frame with 1 row

trial.0 <- try(
               portfolio:::.df.category.mean(data.frame(), by.var = "by.var"),
               silent = TRUE
               )

trial.1 <- try(
               portfolio:::.df.category.mean(data.frame(c("foo"),
                                             by.var = "by.var")),
               silent = TRUE
               )

if(class(trial.0) == "try-error" || class(trial.1) == "try-error"){
  stop("Fails on corner cases!")
}



## Invalid input: list contains non-data.frame objects

trial.2 <- try(
               portfolio:::.df.category.mean(matrix(), array(), data.frame(),
                                             by.var = "by.var"), silent = TRUE
               ) 

if(class(trial.2) == "try-error"){
  stopifnot(
            isTRUE(as.logical(grep("Error.*is\\.data\\.frame",trial.2[1])))
            )
}

## Tests function results against pre-calculated data

stopifnot(
          all.equal(truth$mean.1.2,
                    portfolio:::.df.category.mean(truth$x.1, truth$x.2,
                                                  by.var = "by.var")
                    ),
          all.equal(truth$mean.1.3,
                    portfolio:::.df.category.mean(truth$x.1, truth$x.3,
                                                  by.var = "by.var")
                    ),
          all.equal(truth$mean.1.2.3,
                    portfolio:::.df.category.mean(truth$x.1, truth$x.2, truth$x.3,
                                                  by.var = "by.var")
                    )
          )
