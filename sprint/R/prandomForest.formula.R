##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

## Preprocessing copied shamelessly from serial randomForest implementation
`prandomForest.formula` <-
  function(formula, data = NULL, ..., subset, na.action = na.fail) {
    if (!inherits(formula, "formula"))
      stop("method is only for formula objects")
    m <- match.call(expand.dots = FALSE)
    if (any(c("xtest", "ytest") %in% names(m)))
      stop("xtest/ytest not supported through the formula interface")
    names(m)[2] <- "formula"
    if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    y <- model.response(m)
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    ## Drop any "negative" terms in the formula.
    ## test with:
    ## randomForest(Fertility~.-Catholic+I(Catholic<50),data=swiss,mtry=2)
    m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)),
                     data.frame(m))
    ## if (!is.null(y)) m <- m[, -1, drop=FALSE]
    for (i in seq(along=ncol(m))) {
      if (is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
    }
    ret <- prandomForest(m, y, ...)
    if ( is.null(ret) ) {
      return(ret)
    }
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    ret$call <- cl
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action")))
      ret$na.action <- attr(m, "na.action")
    class(ret) <- c("randomForest.formula", "randomForest")
    return(ret)
  }
