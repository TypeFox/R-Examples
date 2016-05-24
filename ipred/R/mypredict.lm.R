# $Id: mypredict.lm.R,v 1.7 2003/04/02 11:22:49 peters Exp $

mypredict.lm <- function(object, newdata) {

  xn <- as.data.frame(newdata)

  test <- attr(terms(object), "term.labels")
  xn <- xn[,test]

  if (!is.null(nrow(xn))) {
    pred <- rep(NA, nrow(xn))
    names(pred) <- row.names(xn)
  } else {
    pred <- NA
    names(pred) <- "1"
  }

  # evaluate na.omit (delete lines containing NA)

  xnn <- na.omit(xn)

  # attr(xnn, "na.action") returns which na.action is 
  # evaluated, lines and corresponding row.name where NAs occur 

  if(is.null(attr(xnn, "na.action"))) 
    pred <- predict(object, xnn) 
  else 
    pred[-attr(xnn, "na.action")] <- predict(object, xnn)

  pred

}
