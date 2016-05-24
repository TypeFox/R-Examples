ppca.scores.plot <-
function(output, group=FALSE)
{
 q<-output$q
 gpnames<-levels(as.factor(group))
 sigma<-output$sig*solve((t(output$loadings)%*%output$loadings) + (output$sig*diag(ncol(output$loadings))))
 
 if(q==1)
 {
	if((group==FALSE)[1])
    {
        plot(output$scores[,1], rep(0,nrow(output$scores)), yaxt='n', pch=16, font=5, xlab="PC 1", ylab="", col=1)
        legend("topright", bty="n", paste("Variance = ", round(sigma,2), sep=""))
        abline(h=0)
     }else{
         plot(output$scores[,1], rep(0,nrow(output$scores)), yaxt='n', cex=1, pch=rep(16:26)[group], col=group, font=5, xlab="PC 1", ylab="")
         abline(h=0)
         legend("topleft", bty = "n", paste("Group ", c(gpnames), sep=""), col = group, pch = rep(16:26))
         legend("topright", bty="n", paste("Variance = ", round(sigma,2), sep=""))
      }
 }
 
 if(q==2){
     if((group==FALSE)[1])
     {
         plot(output$scores[,1], output$scores[,2], xlim=c(min(output$scores[,1]),max(output$scores[,1])),  ylim=c(min(output$scores[,2]), max(output$scores[,2])+0.5*range(output$scores[,2])[2]), type="n", xlab="PC1", ylab="PC 2")
         for(i in 1:nrow(output$scores))
         {
               points(ellipse(sigma[1:2,1:2], centre=output$scores[, 1:2][i, ], level=0.95), type="l", col="grey50") 
         }
         points(output$scores[,1], output$scores[,2], cex=1, pch=16, font=5, col=1)
         legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = 1)
      }else{
          plot(output$scores[,1], output$scores[,2], xlim=c(min(output$scores[,1]),max(output$scores[,1])), ylim=c(min(output$scores[,2]), max(output$scores[,2])+0.5*range(output$scores[,2])[2]), type="n", xlab="PC1", ylab="PC 2")
          for(i in 1:nrow(output$scores))
          {
               points(ellipse(sigma[1:2,1:2], centre=output$scores[, 1:2][i, ], level=0.95), type="l", col="grey50") 
          }
          points(output$scores[,1], output$scores[,2], cex=1, pch=rep(16:26)[group], col=group, font=5)
          legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          legend("topright", bty = "n", paste("Group ", c(gpnames), sep=""), col = group,lty = rep(0,length(gpnames)),  pch = rep(16:26))
      } # else
 } # if
 
 if(q>2)
 {
 	if((group==FALSE)[1])
    {
    	c1<-1
    	c2<-1
    	for(j in 1:choose(q,(q-2)))
        {
        	c2<-c2+1
        	scores<-cbind(output$scores[,c1], output$scores[, c2])
        	varcov<-matrix(c(sigma[c1,c1], sigma[c1,c2], sigma[c2,c1], sigma[c2,c2]), 2, 2, byrow=TRUE)
           	plot(scores[,1], scores[,2], xlim=c(min(scores[,1])-0.5,max(scores[,1])+0.5), ylim=c(min(scores[,2])-0.5, max(scores[,2])+0.5), type="n", xlab=paste("PC",c1, sep=""), ylab=paste("PC",c2, sep=""))
          	for(l in 1:nrow(scores))
          	{
               points(ellipse(varcov, centre=scores[l, ], level=0.95), type="l", col="grey50") 
          	}
          	points(scores[,1], scores[,2], cex=1, pch=16, font=5, col=1)
          	legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          	if(j < choose(q,(q-2)))
            {
              ask(msg = "Press <RETURN> to view the scores plot for the next pair of dimensions: ")
         	}
         	if(c2 == q){c1<-c1+1; c2<-c1}
        }
     }else{
           c1<-1
    	   c2<-1
    	   for(j in 1:choose(q,(q-2)))
           {
        	 c2<-c2+1
        	 scores<-cbind(output$scores[,c1], output$scores[, c2])
        	 varcov<-matrix(c(sigma[c1,c1], sigma[c1,c2], sigma[c2,c1], sigma[c2,c2]), 2, 2, byrow=TRUE)
           	plot(scores[,1], scores[,2], xlim=c(min(scores[,1])-0.5,max(scores[,1])+0.5), ylim=c(min(scores[,2])-0.5, max(scores[,2])+0.5), type="n", xlab=paste("PC",c1, sep=""), ylab=paste("PC",c2, sep=""))
          	for(l in 1:nrow(scores))
          	{  
               points(ellipse(varcov, centre=scores[l, ], level=0.95), type="l", col="grey50")
          	}
          	points(scores[,1], scores[,2], cex=1, font=5, pch=rep(16:26)[group], col=group)
          	legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          	legend("topright", bty = "n", paste("Group ", c(gpnames), sep=""), col = group,lty = rep(0,length(gpnames)),  pch = rep(16:26))
          	if(j < choose(q,(q-2)))
            {
              ask(msg = "Press <RETURN> to view the scores plot for the next pair of dimensions: ")
         	}
         	if(c2 == q){c1<-c1+1; c2<-c1}
          } #j
        }
 }
} # End of plot.ppca.scores

