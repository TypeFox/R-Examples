plot.tbss <- function(x, first = 2, last = 2, datatype = NULL, main = "The components with most extreme kurtoses", ...){
  if(is.null(datatype)){
    datatype = x$datatype
  }
  
  datatype <- match.arg(datatype, c("iid", "ts"))
  
  switch(datatype,
        iid = {
               pairs(selectComponents(x$S, first = first, last = last), main = main, ...)
               }
        ,
        ts={
               plot.ts(selectComponents(x$S, first = first, last = last), main = main, ...)
             }
        )
}
