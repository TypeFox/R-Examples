plot.MatData <-
function(x,pred.50,pred.95,top.text,...)
    {
     options(warn=-1);
     pred.plot <- min(x$pred):max(x$pred);
     countot   <- vector(mode="numeric",length=length(x$pred)/2);
     prop.obs  <- vector(mode="numeric",length=length(x$pred)/2);
     for (i in 1:length(countot)) 
        {
        countot[i] <- x$count[i*2-1]+x$count[i*2];
        if(countot[i]>0)
           {
            prop.obs[i] <- x$count[i*2]/countot[i];
           }
        else
           {    
            prop.obs[i] <- NA;
           } 
        }
     class(pred.50) <- "numeric";
     class(pred.95) <- "numeric";
     prop.ini <- 1/(1+exp((log(1/19))*((pred.plot-pred.50)/(pred.95-pred.50))));
     plot(pred.plot,prop.obs,main=top.text,ylim=c(0,1.1));
     lines(pred.plot,prop.ini);
    }
