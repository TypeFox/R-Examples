logistic.dma <-
function(x, y, models.which, lambda=0.99, alpha=0.99,autotune=TRUE, 
         initmodelprobs=NULL,initialsamp=NULL) UseMethod("logistic.dma")

