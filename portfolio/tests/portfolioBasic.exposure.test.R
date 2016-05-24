################################################################################
##
## $Id: portfolioBasic.exposure.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests for the portfolioBasic class.
##
################################################################################

library(portfolio)

load("portfolioBasic.exposure.test.RData")

## save(data, exp.1, x, file = "portfolioBasic.exposure.test.RData", compress = TRUE)

## data <- data.frame(id = 1:20, in.var = 1:20)
## data$in.var <- as.numeric(data$in.var)

## ## Construct by.var's.  Note that the assigned vectors are
## ## deliberately recycled.

## data$by.var.1 <- c("1","2")
## data$by.var.2 <- c("1","2","3","4")
## data$by.var.3 <- c(-1,0,0,1)

## x <- new("portfolioBasic", in.var = "in.var", type = "sigmoid", size = 8, data = data)
## x <- create(x)

exp.1.test <- exposure(x, exp.var = c("by.var.1","by.var.2","by.var.3"))

## Function to test equality of data frame components of the exposure
## object's data slot.  Different default handling of row.names in
## 2.4.0 makes it necessary to compare column-by-column.

exp.df.equal <- function(e1, e2) {
  return(all.equal(e1$variable, e2$variable) &&
         all.equal(e1$long,     e2$long) &&
         all.equal(e1$short,    e2$short) &&
         all.equal(e1$exposure, e2$exposure))
}

stopifnot(all.equal(names(exp.1@data), names(exp.1.test@data)),
          exp.df.equal(exp.1@data$numeric,  exp.1.test@data$numeric),
          exp.df.equal(exp.1@data$by.var.1, exp.1.test@data$by.var.1),
          exp.df.equal(exp.1@data$by.var.2, exp.1.test@data$by.var.2))
