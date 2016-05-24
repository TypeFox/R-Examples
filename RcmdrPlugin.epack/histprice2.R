histprice2<- function(inst1,start1="1998-01-01",quot1="Close",end1) {

  library(tseries)
 z <- get.hist.quote(instrument=inst1, start=start1,end=end1,
                     quote=quot1,comp = "m")
y <- as.ts(aggregate(z, as.yearmon, tail, 1))
y.df <- data.frame(y=y,time=time(y))
y.df$x <- ts(y.df[,1])
tsp(y.df$x) <- tsp(y.df[,2])
names(y.df) <- c("data","time","ts")
z.df <- data.frame(ts=y.df$ts)
 return(z.df)
}
