repsaledata <- function(price,timevar,id) {

  data1  <- data.frame(price,timevar,id)
  o <- order(id,timevar,price)
  n = length(o)
  data1 <- data1[o,]
  data0 <- data1[1:(n-1),]
  data1 <- data1[2:n,]
  rsale <- data0$id==data1$id & data0$timevar<data1$timevar
  data0 <- data0[rsale,]
  data1 <- data1[rsale,]
  id <- data1$id
  time0 <- data0$timevar
  time1 <- data1$timevar
  price0 <- data0$price
  price1 <- data1$price
  rdata <- data.frame(id,time0,time1,price0,price1)

  return(rdata)
}


