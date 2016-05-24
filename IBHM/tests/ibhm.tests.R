library(IBHM)
source('utils.R')


(test.Tanh <- function(){
  x <- seq(-6,6,length.out=400)
  y <- tanh(x)
  
  set.seed(123)  
  m <- TrainIBHM(x,y, ConfigureIBHM(scal.candidates='dot.pr', activ.candidates='tanh'))  
  stopifnot(summary(m)$rmse<0.1) 
  cat('tanh - ok\n')
})()

(test.LogSig <- function(){
  x <- seq(-6,6,length.out=400)
  y <- tanh(x)
  
  set.seed(123)  
  m <- TrainIBHM(x,y, ConfigureIBHM(scal.candidates='dot.pr', activ.candidates='logsig'))
  stopifnot(summary(m)$rmse<0.1) 
  cat('logsig - ok\n')
})()


(test.Lin <- function(){
  x <- seq(-6,6,length.out=400)
  y <- 2*x 
  
  set.seed(123)  
  m <- TrainIBHM(x,y, ConfigureIBHM(scal.candidates='dot.pr', activ.candidates='lin'))
  
  stopifnot(summary(m)$rmse<0.1) 
  cat('lin - ok\n')
})()


(test.Tanh2d <- function(){  
  x <- mesh(seq(-3,3,length=50),2)
  y <- tanh(2*(x%*%c(0.1,0.5))+1) + 1
      
  set.seed(123)  
  m <- TrainIBHM(x,y,ConfigureIBHM(scal.candidates='dot.pr', activ.candidates='tanh'))
  
  stopifnot(summary(m)$rmse<0.1) 
  cat('tanh2d - ok\n')
})()

(test.FinalEstimation <- function(){  
  x <- mesh(seq(-3,3,length=50),2)
  y <- tanh(2*(x%*%c(0.1,0.5))+1) + 1
  
  set.seed(123)  
  m <- TrainIBHM(x[1:30,],y[1:30],ConfigureIBHM(scal.candidates='dot.pr', activ.candidates='tanh', final.estimation.x = x, final.estimation.y = y))
  
  stopifnot(summary(m)$rmse<0.1) 
  cat('tanh2d - ok\n')
})()




