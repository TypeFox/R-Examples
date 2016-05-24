dma <-
function(x, y, models.which, lambda=0.99, gamma=0.99, 
         eps=.001/nrow(models.which), delay=1, initialperiod=200) UseMethod("dma")

