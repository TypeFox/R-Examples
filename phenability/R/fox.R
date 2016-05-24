fox <-
function(data, interaction=FALSE){
  a <- nrow(data)
  b <- ncol(data)
  data.m <- as.data.frame(data)
  l.data <- length (data.m)  
  ranks <- matrix(NA,a,b)
  y <- numeric()
  k <- numeric()
  for(i in 1:nrow(data.m)){
    for (j in 1:ncol(data.m)){
      
      ranks[,j] <- rank(-data.m[,j])
    }
    
    y <- which(ranks[i,] <= 3)
    k[i] <- length(y)
    
    
  }
  means <- round(as.numeric(rowMeans(data)),digits=4)
  result <- as.data.frame(cbind(rownames(data),means,k))
  colnames(result) <- c("Gen","Mean", "TOP")
  
  rank.y <- apply(-data,2,rank)
  ranks.sum.y <- apply(rank.y,1,sum)
  sd.rank = round(apply(rank.y,1,sd),digits=4)
  ranks.y = data.frame(rank.y,ranks.sum.y,sd.rank)
  colnames(ranks.y) = c(colnames(rank.y),"Sum", "Sd")
  cor.rank.y <- round(cor(ranks.y, method="pearson"), digits = 4)
  geral.list <- list("Fox"=result,"Ranks"=ranks.y,"Correlations"=cor.rank.y)
  
  if(interaction){
    amb1 = data.frame(data[1:nrow(data),1])
    colnames(amb1) = "amb"
    ambs = amb1
    amb2 <- NULL
    for(j in 2:ncol(data)){
      
      amb2 = data.frame(data[1:nrow(data),j])
      colnames(amb2) = "amb"
      ambs = rbind(ambs, amb2)    
    }
    
    gen <- rep(1:nrow(data),ncol(data))
    env <- rep(1:ncol(data),each=nrow(data))
    intera.data <- data.frame(gen,env,ambs)
    interaction.plot(reorder(factor(intera.data$env),intera.data$amb,mean),
                     intera.data$gen,intera.data$amb,legend = F, type="l",
                     trace.label = deparse(substitute(intera.data$gen)),
                     col = 1:nrow(data),xpd = NULL,xtick = F,cex.axis = 0.6,
                     ylab="Response", 
                     xlab="Environment")
    
  }
  
  top <- as.numeric(result[,3])
  plot(means,top,pch=19,cex=0.5,main="Means x TOP",xlab=expression(Mean[Phenotypic]),
       ylab=expression(TOP[third]),xlim=c(min(means),max(means)),
       ylim=c(min(top),max(top)))
  m <- apply(cbind(means,top),2,mean)
  textxy(means,top,1:a,m=m,cex=1, col="blue")
  origin(m)
  
  return(geral.list)
  
}
