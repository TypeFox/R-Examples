synth.peak <- function(base=0.07,base.time=6, rise.time=5, rise.factor, recession.const=0.2, length.out=240, rez.time=length.out-ceiling(base.time)-ceiling(rise.time)){

   stopifnot(ceiling(base.time)+ceiling(rise.time)+ceiling(rez.time)==ceiling(length.out)) 

   rise.coef <-  log(rise.factor)/(rise.time-1)
   if(base.time <= 0){
      ser.base <- c()
      if(base.time <= -1){
          rise.start <- 1 - base.time
      } else {
          rise.start <- 1
      }
   } else {
      ser.base <- rep(base,ceiling(base.time))
      rise.start <- 1
   }
   if(rise.start > rise.time){
      ser.rise <- c()
   } else {
      ser.rise <- base*exp((rise.start:rise.time-1)*rise.coef)
   }
   ser.rez <- base+(base*exp((rise.time-1)*rise.coef)-base)*exp((1:rez.time)*-recession.const)
   if(rez.time>length.out){
      ser.rez <- ser.rez[-(1:(rez.time-length.out))]
   }
   ans <- c(ser.base,ser.rise,ser.rez)
   stopifnot(length(ans)==ceiling(length.out)) 
   return(ans)
}

