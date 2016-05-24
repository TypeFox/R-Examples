summary.rec.ev.data.sim <-
function(object, ...)
{
   if(!inherits(object, "rec.ev.data.sim")) stop("Wrong data type")
   sub.risk    <- vector()
   num.events  <- vector()
   foltime     <- vector()
   med.foltime <- vector()
   mean.ep.sub <- vector()
   dens.incid  <- vector()
   ep.num      <- vector()
   for (i in 1:attr(object,"ndist"))
   {
     ep.num[i]      <- i
     if (i != 1)
     {
        sub.risk[i]    <- dim(object[object$obs.episode == i-1 & object$status == 1,])[1]
     }else{
        sub.risk[i] <- attr(object,"n")
     }
     num.events[i]  <- as.integer(sum(object$status[object$obs.episode==i]))
     foltime[i]     <- sum(object$time2[object$obs.episode==i])
     med.foltime[i] <- median(object$time2[object$obs.episode==i]) 
     dens.incid[i]  <- num.events[i]/foltime[i]
   }
   ans <- data.frame(ep.num, sub.risk, num.events, 
                     foltime, med.foltime, dens.incid)
  
   ep.num      <- paste (">", attr(object,"ndist"), sep = " ")
   sub.risk     <- dim(object[object$obs.episode == attr(object,"ndist") & object$status == 1,])[1]
   num.events     <- as.integer(sum(object$status[object$obs.episode > attr(object,"ndist")]))
   foltime     <- sum(object$time2[object$obs.episode > attr(object,"ndist")])
   med.foltime     <- median(object$time2[object$obs.episode > attr(object,"ndist")])
   dens.incid     <- num.events/foltime
   
   ans2 <- data.frame(ep.num, sub.risk, num.events, 
                     foltime, med.foltime, dens.incid)
   
   ep.num      <- paste ("All episodes", sep = " ")
   sub.risk     <- attr(object,"n")
   num.events     <- as.integer(sum(object$status))
   foltime     <- sum(object$time2)
   med.foltime     <- median(object$time2)
   dens.incid     <- num.events/foltime
   
   ans3 <- data.frame(ep.num, sub.risk, num.events, 
                      foltime, med.foltime, dens.incid)
   
   ans4 <- rbind(ans,ans2,ans3)
   class(ans4) <- "summary.rec.ev.data.sim"
   return(ans4)
}
