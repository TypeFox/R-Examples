ppcca.scores.plot <-
function(output, Covars, group=FALSE, covarnames=NULL)
{
 q<-output$q
 Covars<-as.matrix(Covars)
 Covars<-standardize(Covars)
 if(colnames(Covars, do.NULL=FALSE)[1] =="col1")
 {
 	colnames(Covars)<-c(paste("Covariate_ ", 1:ncol(Covars), sep=""))
 }
 if(is.null(covarnames) == FALSE)
 {
 	colnames(Covars)<-covarnames
 }
 gpnames<-levels(as.factor(group))
 
 sigma<-output$sig*solve((t(output$loadings)%*%output$loadings) + (output$sig*diag(ncol(output$loadings))))
 
 if(q==1)
  {
	if((group==FALSE)[1])
    {
    	for(i in 1:ncol(Covars))
    	{
         plot(output$scores[,1], rep(0,nrow(output$scores)), yaxt='n', cex=(Covars[,i]+0.5)*2, pch=16, font=5, xlab="PC 1", ylab="", main=colnames(Covars)[i], col=1)
         legend("topright", bty="n", paste("Variance = ", round(sigma,2), sep=""))
         abline(h=0)
         if((ncol(Covars) > 1) & (i < ncol(Covars)))
         {
              ask(msg = "Press <RETURN> to view the scores plot for the next covariate: ")
         }
        }#i
     }else{
     	for(i in 1:ncol(Covars))
     	{
         plot(output$scores[,1], rep(0,nrow(output$scores)), yaxt='n', cex=(Covars[,i]+0.5)*2, pch=16, col=group, font=5, xlab="PC 1", ylab="", main=colnames(Covars)[i])
         legend("topright", bty="n", paste("Variance = ", round(sigma,2), sep=""))
         abline(h=0)
         legend("topleft", bty = "n", paste("Group ", c(gpnames), sep=""), col = group, pch = 16)
         if((ncol(Covars) > 1) & (i < ncol(Covars)))
         {
              ask(msg = "Press <RETURN> to view the scores plot for the next covariate: ")
         }
        }#i
      }
  }
  
  if(q==2)
  {
  	if((group==FALSE)[1])
     {
     	for(i in 1:ncol(Covars))
     	{
         plot(output$scores[,1], output$scores[,2], xlim=c(min(output$scores[,1])-0.5, max(output$scores[,1])+0.5),  ylim=c(min(output$scores[,2])-0.5, max(output$scores[,2])+0.5), xlab="PC1", ylab="PC 2", main=colnames(Covars)[i], cex=(Covars[,i]+0.5)*2, pch=16, font=5, type="n")
         for(j in 1:nrow(output$scores))
         {
               points(ellipse(sigma[1:2,1:2], centre=output$scores[,1:2][j,], level=0.95), type="l", col="grey50")
         }
         points(output$scores[,1], output$scores[,2], cex=(Covars[,i]+0.5)*2, pch=16, font=5, col=1)
         legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
         if((ncol(Covars) > 1) & (i < ncol(Covars)))
         {
              ask(msg = "Press <RETURN> to view the scores plot for the next covariate: ")
         }
        }#i
     }else{
     	for(i in 1:ncol(Covars))
     	{
          plot(output$scores[,1], output$scores[,2], xlim=c(min(output$scores[,1])-0.5, max(output$scores[,1])+0.5),  ylim=c(min(output$scores[,2])-0.5, max(output$scores[,2])+0.5),  xlab="PC1", ylab="PC 2", main=colnames(Covars)[i], cex=(Covars[,i]+0.5)*2, pch=16, font=5, type="n")
         for(j in 1:nrow(output$scores))
         {
               points(ellipse(sigma[1:2,1:2], centre=output$scores[,1:2][j,], level=0.95), type="l", col="grey50")

         }
         points(output$scores[,1], output$scores[,2], cex=(Covars[,i]+0.5)*2, pch=16, font=5, col=group)
         legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          legend("topright", bty = "n", c(paste("Group ", c(gpnames), sep="")), col = group, lty = rep(0,length(gpnames)),  pch = 16)
          if((ncol(Covars) > 1) & (i < ncol(Covars)))
         {
              ask(msg = "Press <RETURN> to view the scores plot for the next covariate: ")
         }
        }#i
       } # else
  } # if

  if(q>2)
  {
  	if((group==FALSE)[1])
    {
        for(i in 1:ncol(Covars))
        {
        	c1<-1
    		c2<-1
    		for(j in 1:choose(q,(q-2)))
        	{
        		c2<-c2+1
        		scores<-cbind(output$scores[,c1], output$scores[, c2])
        		varcov<-matrix(c(sigma[c1,c1], sigma[c1,c2], sigma[c2,c1], sigma[c2,c2]), 2, 2, byrow=TRUE)
           		plot(scores[,1], scores[,2], xlim=c(min(scores[,1])-0.5,max(scores[,1])+0.5), ylim=c(min(scores[,2])-0.5, max(scores[,2])+0.5), type="n", xlab=paste("PC",c1, sep=""), ylab=paste("PC",c2, sep=""), main=colnames(Covars)[i])
          		for(l in 1:nrow(scores))
          		{
               		points(ellipse(varcov, centre=scores[l,], level=0.95), type="l", col="grey50")
          		}
          		points(scores[,1], scores[,2], cex=(Covars[,i]+0.5)*2, pch=16, font=5, col=1)
          		legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          		if(j < choose(q,(q-2)))
            	{
              		ask(msg = "Press <RETURN> to view the scores plot for the next pair of dimensions: ")
         		}
         		if(c2 == q){c1<-c1+1; c2<-c1}
        	}
        	if((ncol(Covars) > 1) & (i < ncol(Covars)))
            {
              ask(msg = "Press <RETURN> to view the scores plot for the next covariate: ")
            }
       } #i
    }else{
       for(i in 1:ncol(Covars))
        {
        	c1<-1
    		c2<-1
    		for(j in 1:choose(q,(q-2)))
        	{
        		c2<-c2+1
        		scores<-cbind(output$scores[,c1], output$scores[, c2])
        		varcov<-matrix(c(sigma[c1,c1], sigma[c1,c2], sigma[c2,c1], sigma[c2,c2]), 2, 2, byrow=TRUE)
           		plot(scores[,1], scores[,2], xlim=c(min(scores[,1])-0.5,max(scores[,1])+0.5), ylim=c(min(scores[,2])-0.5, max(scores[,2])+0.5), type="n", xlab=paste("PC",c1, sep=""), ylab=paste("PC",c2, sep=""), main=colnames(Covars)[i])
          		for(l in 1:nrow(scores))
          		{
               		points(ellipse(varcov, centre=scores[l,], level=0.95), type="l", col="grey50")
          		}
          		points(scores[,1], scores[,2], cex=(Covars[,i]+0.5)*2, pch=16, font=5, col=group)
          		legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          		 legend("topright", bty = "n", paste("Group ", c(gpnames), sep=""), lty = rep(0,length(gpnames)), col = group, pch = 16)
          		if(j < choose(q,(q-2)))
            	{
              		ask(msg = "Press <RETURN> to view the scores plot for the next pair of dimensions: ")
         		}
         		if(c2 == q){c1<-c1+1; c2<-c1}
        	} #j
        	if((ncol(Covars) > 1) & (i < ncol(Covars)))
            {
              ask(msg = "Press <RETURN> to view the scores plot for the next covariate: ")
            }
       } #i
    } #else
  }
} # End of plot.ppcca.scores

