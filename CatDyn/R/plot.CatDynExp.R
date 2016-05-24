plot.CatDynExp <-
function(x,leg.pos,Biom.tstep,Biom.xpos,Biom.ypos,...)
    {
     fleet.name <- x$Properties$Fleets$Fleet
     period     <- x$Model$Results[,1]
     tstep      <- x$Properties$Units[1]
     catunits   <- paste("Catch (",x$Properties$Units[4],")",sep="")
     index      <- 1:length(x$Model$Results[,dim(x$Model$Results)[2]])
     Biom       <- if(Biom.tstep == "fin")
                     {
                      round(mean(tail(x$Model$Results[,"Pred. Biomass (tonnes)"],5)))
                     }
                   else if(Biom.tstep=="mid")
                     {
                      round(mean(x$Model$Results[,"Pred. Biomass (tonnes)"][c(floor(median(index-2)),floor(median(index-1)),floor(median(index)),floor(median(index+1)),floor(median(index+2)))]))
                     }
                   else
                     {
                      round(mean(head(x$Model$Results[,"Pred. Biomass (tonnes)"],5)))
                     }
     options(warn=-1)
     for(i in 1:length(fleet.name))
       {
       if(any(is.nan(x$Model$Results[,i*4]))) 
         {stop(paste("Change initial parameter values to increase predicted catch of ",fleet.name[i], " fleet",sep=""))}
       obscat <- x$Model$Results[,i*3+3*(i-1)]
       modcat <- x$Model$Results[,i*4+2*(i-1)]
       resids <- x$Model$Results[,i*7-1*(i-1)]
       par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(4,4,2,2))
       plot(x=period,y=obscat,pch=19, xlab=gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", tstep, perl=TRUE),ylab=catunits,ylim=c(0,max(obscat,modcat)),main="",cex=0.75)
       lines(x=period,y=modcat,col="black",lwd=2)
       legend(leg.pos,c("Observed Catch", "Predicted Catch"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
       if(x$Properties$Units[3] != "ind")
         {
          text(x=Biom.xpos*max(period), y=Biom.ypos*max(obscat), adj=0,lab=list(bquote("Biom (tons)" ==.(Biom))))
         }
       if(x$Model$Type[i] != 0)
         {
          if(i == 1)
            {
             if(x$Model$Type[i] > 0)
               {
                points(x=x$Model$Dates[2:(x$Model$Type[i]+1)],
                       y=rep(0,x$Model$Type[i]),
                       pch=10,cex=3)
               }
             else
               {
                points(x=x$Model$Dates[seq(2,(2*abs(x$Model$Type[i])+1),2)],
                       y=rep(0,abs(x$Model$Type[i])),
                       pch=10,cex=1,col="red")
                points(x=x$Model$Dates[seq(3,(2*abs(x$Model$Type[i])+1),2)],
                       y=rep(0,abs(x$Model$Type[i])),
                       pch=10,cex=1,col="blue")
               }
            }
          else
            {
             points(x=x$Model$Dates[(2+x$Model$Type[i-1]):(1+x$Model$Type[i-1]+x$Model$Type[i])],
                    y=rep(0,x$Model$Type[i]),
                    pch=10,cex=3)
            }
         }
       hist(x=resids,main="",xlab="Deviance Residuals",ylab="Frequency")
       plot(x=period,y=resids,xlab=gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", tstep, perl=TRUE), ylab="Deviance Residuals",pch=1,type="n")
       abline(h=0,lwd=2)
       text(x=period,y=resids,lab=format(period),cex=0.75)
       qqplot(x=obscat,y=modcat,xlab=paste("Observed Catch (",x$Properties$Units[4],")", sep=""),ylab=paste("Predicted Catch (",x$Properties$Units[4],")", sep=""),pch=1)
       abline(a=0,b=1,lwd=2)
       mtext(side=3,outer=TRUE,text=paste(fleet.name[i]," - ",x$Model$Type[i]," - ",gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", x$Model$Distribution[i], perl=TRUE),sep=""))
       options(warn=0)
       devAskNewPage(ask=TRUE)
       }
     devAskNewPage(ask=FALSE)
    }
