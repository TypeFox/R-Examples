# File cec2013.R
# Part of the cec2013 R package, http://www.rforge.net/cec2013/ ; 
#                                http://cran.r-project.org/web/packages/cec2013
# Copyright 2013 Mauricio Zambrano-Bigiarini & Yasser Gonzalez-Fernandez
# Distributed under GPL 3 or later

################################################################################
#                                 CEC 2013                                     #
################################################################################
# Purpose   : Evaluate a CEC-2013 benchamark function on a user-defined para-  #
#             meter set                                                        #
################################################################################
# i: integer in [1, 28], with the number of the CEC2013 benchmark function to  #
#    be evalauated                                                             #
# x: numeric, with the parameter set to be evaluated in the benchmark function #
#    Its length MUST be in [2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]     #
################################################################################

cec2013 <- function (i, x) {

  if (missing(i)) stop("Missing argument; 'i' has to be provided !")

  if (missing(x)) stop("Missing argument; 'x' has to be provided !")

  if (is.numeric(i) && i >= 1 && i <= 28) {
    if (is.vector(x)) {
        row <- 1; col <- length(x)
    } else if (is.matrix(x)) {
         row <- nrow(x); col <- ncol(x)
      } else {
          stop("x should be a vector or a matrix")
        } # ELSE end

        if (!(col %in% c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))) {
            stop("Invalid argument: Only 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90 and 100 dimensions/variables are allowed !")
        }
        extdatadir <- system.file("extdata", package = "cec2013")
        f <- .C("cec2013", extdatadir = as.character(extdatadir), 
                i = as.integer(i), x = as.double(x), row = as.integer(row),
                col = as.integer(col), f = double(row), 
                PACKAGE = "cec2013")$f
  } else stop("Invalid argument: 'i' should be an integer between 1 and 28 !")
  
  return(f)
}
