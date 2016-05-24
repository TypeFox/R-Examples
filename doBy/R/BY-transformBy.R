

transformBy <- function (formula, data, ...) {

  transform2<- function (data, ...) 
    {
      e <- eval(substitute(list(...)), data, parent.frame())
      tags <- names(e)
      inx <- match(tags, names(data))
      matched <- !is.na(inx)
      if (any(matched)) {
        data[inx[matched]] <- e[matched]
        data <- data.frame(data)
      }
      if (!all(matched)) {
        for (i in 1:length(e))
          data[,names(e)[i]] <- e[i]
      }
      return(data)
    }

  ddd <- splitBy(formula, data=data, drop=TRUE)
  ee <- lapply(ddd, function(d){
    transform2(d, ...)
  })
  do.call("rbind",ee)
}
