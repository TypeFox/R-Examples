################################################################################
#  Prediction of frailties                                                     #
################################################################################
#                                                                              #
#  Computes the prediction of the fraities as                                  #
#                                                                              #
#  Its only parameter is                                                       #
#   - object  : the fitted model, object of class 'parfm'                      #
#                                                                              #
#   Date: February 02, 2012                                                    #
#   Last modification on: May 15, 2012                                         #
################################################################################

predict.parfm <- function(object,
                          ...) {
# Frailty distribution
  if (attributes(object)$frailty == "none")
    stop("The model 'object' is is a simple Cox modelm with no frailties!")
  frailty <- eval(parse(text=paste("fr", attributes(object)$frailty, sep=".")))
  frPar  <- c(rownames(object)[1],
              object[1, 1])
  
  # Baseline hazard
  dist <- eval(parse(text=attributes(object)$dist))
  
# Data needed for the derivatives of the  Laplace transform 
  cumhaz <- attributes(object)$cumhaz
  di <- attributes(object)$di
  clusters <- names(attr(object, "di"))

  res <- sapply(clusters, FUN=function(h) {
    exp(diff(sapply(
      paste("frailty(k=attributes(object)$di['", h, "']+", 0:1,
                    ", s=cumhaz['", h, "'], ",
            paste(frPar, collapse="="), ", ",
            ifelse(attributes(object)$frailty == "possta", 
                   paste("Omega=Omega(D=max(di)+1, ",
                         "correct=", attr(object, 'correct'), 
                         ", nu=", frPar[2] ,")",
                         ", correct=", attr(object, 'correct'), ", ",
                         sep=""),
                   ""),
            "what='logLT'",
            ")", sep=""),
      function(x) eval(parse(text=x)),
   USE.NAMES=FALSE)))
  }, USE.NAMES=FALSE)
  
  class(res) <- "predict.parfm"
  attr(res, "clustname") <-attr(object, "clustname")
  attr(res, "frailty") <-attr(object, "frailty")
  attr(res, "dist") <-attr(object, "dist")
  return(res)
}
