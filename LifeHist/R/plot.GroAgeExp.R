plot.GroAgeExp <-
function(x, ...)
    {
     spec.name  <- x$Properties$Species
     sex        <- x$Properties$Sex
     model      <- x$Properties$Model
     distr      <- x$Properties$Distr
     options(warn=-1)
     if(x$Properties$Sex == "Males" | x$Properties$Sex == "Females" | x$Properties$Sex == "Pooled")
       {
        y      <- x$Results[ do.call(order, list(x$Results$age)), ]
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
     else if(x$Properties$Sex == "Both" | x$Properties$Sex == "Total")
       {
        y          <- vector("list",2)
        y$Females  <- x$Females[ do.call(order, list(x$Females$age)), ]
        y$Males    <- x$Males[ do.call(order, list(x$Males$age)), ]
        age.fem    <- y$Females$age
        obslen.fem <- y$Females$obslen
        modlen.fem <- y$Females$predlen
        resids.fem <- y$Females$resid
        age.mal    <- y$Males$age
        obslen.mal <- y$Males$obslen
        modlen.mal <- y$Males$predlen
        resids.mal <- y$Males$resid
        #layout(matrix(c(1,1,2,2,3,4,5,6),2,4, byrow = TRUE))
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
        mtext(side=3,outer=TRUE,text=paste(spec.name," / ", model," / ", distr, sep=""))
       }
    }
