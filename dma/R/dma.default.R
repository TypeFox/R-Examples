dma.default <-
function(x, y, models.which, lambda=0.99, gamma=0.99, 
 eps=.001/nrow(models.which), delay=1, initialperiod=200) {
    
    #run 
    est<-makf4(x, y, models.which, lambda, gamma,eps, delay, initialperiod)
    
    #define outputs
    est$fitted.values<-est$yhat.ma
    est$residuals<-y-est$yhat.ma
    
    #define as class to allow other commands later
    class(est)<-"dma"
    
    est
}

