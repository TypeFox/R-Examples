"add.log.axis" <-
function(make.labs=FALSE, logs=c(2, 3, 4, 5, 6, 7, 8, 9), side=1,
         two.sided=FALSE, label=NULL, ...) {

   ifelse(side == 2 | side == 4, interval <- par()$usr[3:4],
                                 interval <- par()$usr[1:2])

   interval <- 10^interval

   other.side <- switch(as.character(side), "1"=3, "2"=4, "3"=1, "4"=1)
   interval <- sort(interval)
   min <- interval[1];    max <- interval[2]
   log.beg  <- as.integer(log10(min)) - 1
   log.end  <- as.integer(log10(max)) + 1
   the.logs <- vector(mode="numeric")
   for(cycle in log.beg:log.end) {
      the.logs.in.cycle <- logs*10^(cycle)
      the.logs <- c(the.logs, the.logs.in.cycle)
   }
   if(make.labs) {
      Axis(at=the.logs, labels=the.logs, side=side, tcl=0, ...)
      mtext(label, line=2, side=side)
   } else {
      Axis(at=the.logs, labels=NA, side=side, ...)
      if(two.sided) Axis(at=the.logs,  labels=NA, side=other.side, ...)
   }
}

