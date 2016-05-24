plot.Ord.logreg <-
function(x, one.plot=TRUE,...)
{
 if (class(x)!="Ord.logreg")
stop("x is not of class Ord.logreg")
 model<-x$model
 cat<-length(model)
 if (one.plot==TRUE)
  {
  if(cat<4)
   {
   par(mfrow=c(1,3))
   for (i in 1:cat)
     {
     plot(model[[i]], title=FALSE)
     mtext(paste("Category", i, "Model"), side=3, line=1)
     }
   }
  if(cat==4)
   {
   par(mfrow=c(2,2))
   for (i in 1:cat)
     {
     plot(model[[i]], title=FALSE)
     mtext(paste("Category", i, "Model"), side=3, line=1)
     }
   }
  if(cat>4)
   {
   par(mfrow=c(2,2), ask=TRUE)
   for (i in 1:4)
     {
     plot(model[[i]], title=FALSE)
     mtext(paste("Category", i, "Model"), side=3, line=1)
     }
    plot.new()
    par(mfrow=c(2,2))
    for(j in 5:cat)
     {
     plot(model[[j]], title=FALSE)
     mtext(paste("Category", j, "Model"), side=3, line=1)
     } 
   }
 }
 if (one.plot==FALSE)
  {
  par(ask=TRUE)
  for (i in 1:cat)
     {
     plot.new()
     plot(model[[i]], title=FALSE)
     mtext(paste("Category", i, "Model"), side=3, line=1)
     }
  }
}
