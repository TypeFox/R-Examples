plot.GroAgeMod <-
function(x, ...)
    {
     spec.name  <- x$Properties$Species
     sex        <- x$Properties$Sex
     model      <- x$Model$Growth$Type
     distr      <- x$Model$Distribution$Distr
     method     <- x$Model$Method
     options(warn=-1)
     if(sex == "Males" | sex == "Females" | sex == "Pooled")
       {
        y      <- x$Predictions[ do.call(order, list(x$Predictions$age)), ]
        age    <- y$age
        obslen <- y$obslen
        modlen <- y$predlen
        resids <- y$resid
        par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(4,4,2,2))
        plot(x=age,y=obslen,pch=19,
             xlab=paste("Age (", x$Properties$Units.Age, ")", sep=""),
             ylab=paste(x$Properties$Length.Type, " Length (", x$Properties$Units.Length, ")", sep=""),
             ylim=c(0,max(obslen,modlen)), main="", cex=0.75)
        lines(x=age,y=modlen,col="black",lwd=2)
        legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
        #
        hist(x=resids,main="",xlab="Deviance Residuals",ylab="Frequency")
        #
        plot(x=age,y=resids,xlab=paste("Age (", x$Properties$Units.Age, ")", sep=""),
             ylab="Deviance Residuals",pch=1)
        abline(h=0,lwd=2)
        #text(x=age,y=resids,lab=format(age),cex=0.75)
        #
        qqplot(x=obslen,y=modlen,
               xlab=paste("Observed ",x$Properties$Length.Type, " Length (", x$Properties$Units.Length, ")", sep=""),
               ylab=paste("Predicted ",x$Properties$Length.Type, " Length (", x$Properties$Units.Length, ")", sep=""),
               pch=1)
        abline(a=0,b=1,lwd=2)
        mtext(side=3,outer=TRUE,text=paste(spec.name," - ", sex, " - ", model," - ", distr, sep=""))
       }
     else if(sex == "Both" | sex == "Total")
       {
        if(length(x$Predictions) == 3)
          {
           y          <- vector("list",3)
           y$Females  <- x$Predictions$Females[ do.call(order, list(x$Predictions$Females$age)), ]
           y$Males    <- x$Predictions$Males[ do.call(order, list(x$Predictions$Males$age)), ]
           y$Pooled   <- x$Predictions$Pooled[ do.call(order, list(x$Predictions$Pooled$age)), ]
           age.fem    <- y$Females$age
           obslen.fem <- y$Females$obslen
           modlen.fem <- y$Females$predlen
           resids.fem <- y$Females$resid
           age.mal    <- y$Males$age
           obslen.mal <- y$Males$obslen
           modlen.mal <- y$Males$predlen
           resids.mal <- y$Males$resid
           age.pld    <- y$Pooled$age
           obslen.pld <- y$Pooled$obslen
           modlen.pld <- y$Pooled$predlen
           resids.pld <- y$Pooled$resid
           #
           par(mfrow=c(3,3),oma=c(2,2,2,1), mar=c(4,4,2,2))
           plot(x=age.fem,y=obslen.fem,pch=19,
                xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                ylab=paste(x$Properties$Length.Type, " Length Females (", x$Properties$Units.Length, ")", sep=""),
                ylim=c(0,max(obslen.fem,modlen.fem)), main="", cex=0.75)
           lines(x=age.fem,y=modlen.fem,col="black",lwd=2)
           legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
           #
           plot(x=age.mal,y=obslen.mal,pch=19,
                xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                ylab=paste(x$Properties$Length.Type, " Length Males (", x$Properties$Units.Length, ")", sep=""),
                ylim=c(0,max(obslen.mal,modlen.mal)), main="", cex=0.75)
           lines(x=age.mal,y=modlen.mal,col="black",lwd=2)
           legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
           plot(x=age.pld,y=obslen.pld,pch=19,
                xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                ylab=paste(x$Properties$Length.Type, " Length Pooled (", x$Properties$Units.Length, ")", sep=""),
                ylim=c(0,max(obslen.pld,modlen.pld)), main="", cex=0.75)
           lines(x=age.pld,y=modlen.pld,col="black",lwd=2)
           legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
           #
           hist(x=resids.fem,main="",xlab="Deviance Residuals",ylab="Frequency Females")
           #
           hist(x=resids.mal,main="",xlab="Deviance Residuals",ylab="Frequency Males")
           #
           hist(x=resids.pld,main="",xlab="Deviance Residuals",ylab="Frequency Pooled")
           #
           plot(x=age.fem,y=resids.fem,xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                ylab="Deviance Residuals",pch=1)
           abline(h=0,lwd=2)
           #text(x=age,y=resids,lab=format(age),cex=0.75)
           #
           plot(x=age.mal,y=resids.mal,xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                ylab="Deviance Residuals",pch=1)
           abline(h=0,lwd=2)
           #
           plot(x=age.pld,y=resids.pld,xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                ylab="Deviance Residuals",pch=1)
           abline(h=0,lwd=2)
           #
           mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, " / ", method, sep=""))
          }
        else if(length(x$Predictions) == 2)
          {
           if(names(x$Predictions) == c("Females", "Males"))
             {
              y          <- vector("list",2)
              y$Females  <- x$Predictions$Females[ do.call(order, list(x$Predictions$Females$age)), ]
              y$Males    <- x$Predictions$Males[ do.call(order, list(x$Predictions$Males$age)), ]
              age.fem    <- y$Females$age
              obslen.fem <- y$Females$obslen
              modlen.fem <- y$Females$predlen
              resids.fem <- y$Females$resid
              age.mal    <- y$Males$age
              obslen.mal <- y$Males$obslen
              modlen.mal <- y$Males$predlen
              resids.mal <- y$Males$resid
              par(mfrow=c(3,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              plot(x=age.fem,y=obslen.fem,pch=19,
                   xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Females (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.fem,modlen.fem)), main="", cex=0.75)
              lines(x=age.fem,y=modlen.fem,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              plot(x=age.mal,y=obslen.mal,pch=19,
                   xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Males (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.mal,modlen.mal)), main="", cex=0.75)
              lines(x=age.mal,y=modlen.mal,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              hist(x=resids.fem,main="",xlab="Deviance Residuals",ylab="Frequency Females")
              #
              hist(x=resids.mal,main="",xlab="Deviance Residuals",ylab="Frequency Males")
              #
              plot(x=age.fem,y=resids.fem,xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #text(x=age,y=resids,lab=format(age),cex=0.75)
              #
              plot(x=age.mal,y=resids.mal,xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, " / ", method, sep=""))
             }
           else if(names(x$Predictions) == c("Females", "Pooled"))
             {
              y          <- vector("list",2)
              y$Females  <- x$Predictions$Females[ do.call(order, list(x$Predictions$Females$age)), ]
              y$Pooled   <- x$Predictions$Pooled[ do.call(order, list(x$Predictions$Pooled$age)), ]
              age.fem    <- y$Females$age
              obslen.fem <- y$Females$obslen
              modlen.fem <- y$Females$predlen
              resids.fem <- y$Females$resid
              age.pld    <- y$Pooled$age
              obslen.pld <- y$Pooled$obslen
              modlen.pld <- y$Pooled$predlen
              resids.pld <- y$Pooled$resid
              par(mfrow=c(3,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              plot(x=age.fem,y=obslen.fem,pch=19,
                   xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Females (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.fem,modlen.fem)), main="", cex=0.75)
              lines(x=age.fem,y=modlen.fem,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              plot(x=age.pld,y=obslen.pld,pch=19,
                   xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Pooled (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.pld,modlen.pld)), main="", cex=0.75)
              lines(x=age.pld,y=modlen.pld,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              hist(x=resids.fem,main="",xlab="Deviance Residuals",ylab="Frequency Females")
              #
              hist(x=resids.pld,main="",xlab="Deviance Residuals",ylab="Frequency Pooled")
              #
              plot(x=age.fem,y=resids.fem,xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              plot(x=age.pld,y=resids.pld,xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, " / ", method, sep=""))
             }
           else if(names(x$Predictions) == c("Males", "Pooled"))
             {
              y          <- vector("list",2)
              y$Males    <- x$Predictions$Males[ do.call(order, list(x$Predictions$Males$age)), ]
              y$Pooled   <- x$Predictions$Pooled[ do.call(order, list(x$Predictions$Pooled$age)), ]
              age.mal    <- y$Males$age
              obslen.mal <- y$Males$obslen
              modlen.mal <- y$Males$predlen
              resids.mal <- y$Males$resid
              age.pld    <- y$Pooled$age
              obslen.pld <- y$Pooled$obslen
              modlen.pld <- y$Pooled$predlen
              resids.pld <- y$Pooled$resid
              par(mfrow=c(3,2),oma=c(2,2,2,1), mar=c(4,4,2,2))
              #
              plot(x=age.mal,y=obslen.mal,pch=19,
                   xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Males (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.mal,modlen.mal)), main="", cex=0.75)
              lines(x=age.mal,y=modlen.mal,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              plot(x=age.pld,y=obslen.pld,pch=19,
                   xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Pooled (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.pld,modlen.pld)), main="", cex=0.75)
              lines(x=age.pld,y=modlen.pld,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              hist(x=resids.mal,main="",xlab="Deviance Residuals",ylab="Frequency Males")
              #
              hist(x=resids.pld,main="",xlab="Deviance Residuals",ylab="Frequency Pooled")
              #
              plot(x=age.mal,y=resids.mal,xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              plot(x=age.pld,y=resids.pld,xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, " / ", method, sep=""))
             }
          }
        else if(length(x$Predictions) == 1)
          {
           if(names(x$Predictions) == "Females")
             {
              y          <- vector("list",1)
              y$Females  <- x$Predictions$Females[ do.call(order, list(x$Predictions$Females$age)), ]
              age.fem    <- y$Females$age
              obslen.fem <- y$Females$obslen
              modlen.fem <- y$Females$predlen
              resids.fem <- y$Females$resid
              par(mfrow=c(3,1),oma=c(2,2,2,1), mar=c(4,4,2,2))
              plot(x=age.fem,y=obslen.fem,pch=19,
                   xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Females (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.fem,modlen.fem)), main="", cex=0.75)
              lines(x=age.fem,y=modlen.fem,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              hist(x=resids.fem,main="",xlab="Deviance Residuals",ylab="Frequency Females")
              #
              plot(x=age.fem,y=resids.fem,xlab=paste("Age Females (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, " / ", method, sep=""))
             }
           else if(names(x$Predictions) == "Males")
             {
              y          <- vector("list",1)
              y$Males    <- x$Predictions$Males[ do.call(order, list(x$Predictions$Males$age)), ]
              age.mal    <- y$Males$age
              obslen.mal <- y$Males$obslen
              modlen.mal <- y$Males$predlen
              resids.mal <- y$Males$resid
              par(mfrow=c(3,1),oma=c(2,2,2,1), mar=c(4,4,2,2))
              #
              plot(x=age.mal,y=obslen.mal,pch=19,
                   xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Males (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.mal,modlen.mal)), main="", cex=0.75)
              lines(x=age.mal,y=modlen.mal,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              hist(x=resids.mal,main="",xlab="Deviance Residuals",ylab="Frequency Males")
              #
              plot(x=age.mal,y=resids.mal,xlab=paste("Age Males (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, " / ", method, sep=""))
             }
           else if(names(x$Predictions) == "Pooled")
             {
              y          <- vector("list",1)
              y$Pooled   <- x$Predictions$Pooled[ do.call(order, list(x$Predictions$Pooled$age)), ]
              age.pld    <- y$Pooled$age
              obslen.pld <- y$Pooled$obslen
              modlen.pld <- y$Pooled$predlen
              resids.pld <- y$Pooled$resid
              par(mfrow=c(3,1),oma=c(2,2,2,1), mar=c(4,4,2,2))
              plot(x=age.pld,y=obslen.pld,pch=19,
                   xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                   ylab=paste(x$Properties$Length.Type, " Length Pooled (", x$Properties$Units.Length, ")", sep=""),
                   ylim=c(0,max(obslen.pld,modlen.pld)), main="", cex=0.75)
              lines(x=age.pld,y=modlen.pld,col="black",lwd=2)
              legend("bottomright",c("Observed Length", "Predicted Length"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
              #
              hist(x=resids.pld,main="",xlab="Deviance Residuals",ylab="Frequency Pooled")
              #
              plot(x=age.pld,y=resids.pld,xlab=paste("Age Pooled (", x$Properties$Units.Age, ")", sep=""),
                   ylab="Deviance Residuals",pch=1)
              abline(h=0,lwd=2)
              #
              mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, " / ", method, sep=""))
             }
          }
       }
    }
