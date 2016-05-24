mppca.scores.plot <-
function(output, group=FALSE, gplegend=TRUE)
{
  q<-output$q
  g<-output$g
  scores<-output$scores
  sig<-output$sig
  loadings<-output$loadings
  gpnames<-levels(as.factor(group))
  
   if(q==1)
   {
	 if((group==FALSE)[1])
     {
      	  for(k in 1:g)
  		  {
              plot(scores[[k]][,1], rep(0,nrow(scores[[k]])), yaxt='n', cex=1.3, pch=16, col=1, font=5, xlab="PC 1", ylab="", main=paste("Group ", k, sep=""))
              sigma<-(sig)*solve((t(loadings[,,k])%*%loadings[,,k]) + (sig*diag(q)))
              legend("topright", bty="n", paste("Variance = ", round(sigma,2), sep=""))
              abline(h=0)
              if((g > 1) & (k < g))
              {
               ask(msg = "Press <RETURN> to view the scores plot for the next group: ")
              }
          } #k
      }#if
      else{
      	for(k in 1:g)
  		{
              plot(scores[[k]][,1], rep(0,nrow(scores[[k]])), yaxt='n', cex=1.3, pch=rep(15:26)[group[output$clustering==k]], col=group[output$clustering==k], font=5, xlab="PC 1", ylab="", main=paste("Group ", k, sep=""))
              sigma<-(sig)*solve((t(loadings[,,k])%*%loadings[,,k]) + (sig*diag(q)))
              legend("topright", bty="n", paste("Variance = ", round(sigma,2), sep=""))
              abline(h=0)
              if(gplegend == TRUE)
              {
               legend("topleft", bty = "n", paste("Treatment group ", c(gpnames), sep=""), col = 1:g, pch = rep(15:26)[as.numeric(names(table(group)))])
              }
             if((g > 1) & (k < g))
             {
              ask(msg = "Press <RETURN> to view the scores plot for the next group: ")
             }
          } #k
          }#else
     }#if

     if(q==2)
     {
         if((group==FALSE)[1])
         {
         	for(k in 1:g)
         	{
              plot(scores[[k]][,1], scores[[k]][,2], xlim=c(min(scores[[k]][,1])-1,max(scores[[k]][,1])+1),  ylim=c(min(scores[[k]][,2])-1, 
              max(scores[[k]][,2])+1), type="n", xlab="PC1", ylab="PC 2", main=paste("Group ", k, sep=""))
              sigma<-(sig)*solve((t(loadings[,,k])%*%loadings[,,k]) + (sig*diag(q)))
              for(i in 1:nrow(scores[[k]]))
              { 
                  points(ellipse(sigma[1:2,1:2], centre=scores[[k]][,1:2][i,], level=0.95), type="l", col="grey50")
              }#i
              points(scores[[k]][,1], scores[[k]][,2], cex=1.3, pch=16, font=5, col=1)   
              legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
             if((g > 1) & (k < g))
             {
              ask(msg = "Press <RETURN> to view the scores plot for the next group: ")
             }
          } #k
         }else{
         	for(k in 1:g)
         	{
                  plot(scores[[k]][,1], scores[[k]][,2], xlim=c(min(scores[[k]][,1])-1,max(scores[[k]][,1])+1), ylim=c(min(scores[[k]][,2])-1, 
                  max(scores[[k]][,2])+1), type="n", xlab="PC1", ylab="PC 2", main=paste("Group ", k, sep=""))
                  sigma<-(sig)*solve((t(loadings[,,k])%*%loadings[,,k]) + (sig*diag(q)))
                  for(i in 1:nrow(scores[[k]]))
                  { 
                     points(ellipse(sigma[1:2,1:2], centre=scores[[k]][,1:2][i,], level=0.95), type="l", col="grey50")
                  }#i
                  points(scores[[k]][,1], scores[[k]][,2], cex=1, pch=rep(15:26)[group[output$clustering==k]], col=group[output$clustering==k], font=5)
                  legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
             if(gplegend == TRUE)
             {
             	legend("topright", bty = "n", c(paste("Treatment Group ", c(gpnames), sep="")), col = 1:g, 
                  lty = c(rep(0,length(gpnames))),  pch = rep(15:26)[as.numeric(names(table(group)))])
             }
             if((g > 1) & (k < g))
             {
              ask(msg = "Press <RETURN> to view the scores plot for the next group: ")
             }
          } #k
         } # else
     } # if
     
     if(q>2)
     {        
          if((group==FALSE)[1])
          {
           for(k in 1:g)
           {
           	sigma<-(sig)*solve((t(loadings[,,k])%*%loadings[,,k]) + (sig*diag(q)))
            c1<-1
    	    c2<-1
    	    for(j in 1:choose(q,(q-2)))
            {
        	  c2<-c2+1
        	  tempscores<-cbind(scores[[k]][,c1], scores[[k]][, c2])
        	  varcov<-matrix(c(sigma[c1,c1], sigma[c1,c2], sigma[c2,c1], sigma[c2,c2]), 2, 2, byrow=TRUE)
           	plot(tempscores[,1], tempscores[,2], xlim=c(min(tempscores[,1])-1,max(tempscores[,1])+1), ylim=c(min(tempscores[,2])-1, max(tempscores[,2])+1), type="n", xlab=paste("PC",c1, sep=""), ylab=paste("PC",c2, sep=""), main=paste("Group ", k, sep=""))
          	for(l in 1:nrow(tempscores))
          	{
               points(ellipse(varcov, centre=tempscores[l,], level=0.95), type="l", col="grey50")
          	}
          	points(tempscores[,1], tempscores[,2], cex=1, pch=16, font=5, col=1)
          	legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          	if(j < choose(q,(q-2)))
            {
              ask(msg = "Press <RETURN> to view the scores plot for the next pair of dimensions: ")
         	}
         	if(c2 == q){c1<-c1+1; c2<-c1}
           } #j
           if((g > 1) & (k < g))
           {
              ask(msg = "Press <RETURN> to view the scores plot for the next group: ")
            }
          } #k
         }else{
         	for(k in 1:g)
         	{
            sigma<-(sig)*solve((t(loadings[,,k])%*%loadings[,,k]) + (sig*diag(q)))
            c1<-1
    	    c2<-1
    	    for(j in 1:choose(q,(q-2)))
            {
        	  c2<-c2+1
        	  tempscores<-cbind(scores[[k]][,c1], scores[[k]][, c2])
        	  varcov<-matrix(c(sigma[c1,c1], sigma[c1,c2], sigma[c2,c1], sigma[c2,c2]), 2, 2, byrow=TRUE)
           	plot(tempscores[,1], tempscores[,2], xlim=c(min(tempscores[,1])-1,max(tempscores[,1])+1), ylim=c(min(tempscores[,2])-1, max(tempscores[,2])+1), type="n", xlab=paste("PC",c1, sep=""), ylab=paste("PC",c2, sep=""), main=paste("Group ", k, sep=""))
          	for(l in 1:nrow(tempscores))
          	{
          		points(ellipse(varcov, centre=tempscores[l,], level=0.95), type="l", col="grey50")
          	}
          	points(tempscores[,1], tempscores[,2], cex=1, pch=rep(15:26)[group[output$clustering==k]], font=5, col=group[output$clustering==k])
          	legend("topleft", bty = "n", c("95% Posterior Set"), col = "grey50", lty = c(1))
          	if(gplegend == TRUE)
            {
             	legend("topright", bty = "n", c(paste("Treatment Group ", c(gpnames), sep="")), col = 1:g, lty = c(rep(0,length(gpnames))),  pch = rep(15:26)[as.numeric(names(table(group)))])
            }
          	if(j < choose(q,(q-2)))
            {
              ask(msg = "Press <RETURN> to view the scores plot for the next pair of dimensions: ")
         	}
         	if(c2 == q){c1<-c1+1; c2<-c1}
           } #j
           if((g > 1) & (k < g))
           {
              ask(msg = "Press <RETURN> to view the scores plot for the next group: ")
            }
          } #k
        }#else
     }#if
  
} # End of plot.mppca.scores

