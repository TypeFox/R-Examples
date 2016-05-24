kang <-
function (data,interaction=FALSE){
  a <- nrow(data)
  b <- ncol(data)
  data.m <- as.matrix(data)
  l.data <- length (data.m)
  yield <- numeric()
  vars <- numeric()
  
  for (i in 1:nrow(data.m)){
    
    yield[i] <- sum(data.m[i,])  
    vars[i] <- var(data.m[i,])
  }
  
  rank.y <- rank(-yield)
  rank.v <- rank(vars)
  rank.sum <- rank.y + rank.v
  rank.sum <- as.data.frame(rank.sum)
  
  means <- round(as.vector(rowMeans(data)),digits=4)
  result <- as.data.frame(cbind(rownames(data),means,rank.y,rank.v,rank.sum))
  colnames(result) <- c("Gen","Mean","Rank.Y","Rank.VAR","rank.sum")
  
  rank.y <- apply(-data,2,rank)
  ranks.sum.y <- apply(rank.y,1,sum)
  sd.rank = round(apply(rank.y,1,sd),digits=4)
  ranks.y = data.frame(rank.y,ranks.sum.y,sd.rank)
  colnames(ranks.y) = c(colnames(rank.y),"Sum", "Sd")
  cor.rank.y <- round(cor(ranks.y, method="pearson"), digits = 4)
  geral.list <- list("Kang"=result,"Ranks"=ranks.y,"Correlations"=cor.rank.y)
  
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
  
  rs <- as.numeric(result[,5])
  plot(result[,2],result[,5],pch=19,cex=0.5,main="Means x Rank Sum",
       xlab=expression(Mean[phenotypic]),ylab=expression(rs["(rank sum)"]),
       xlim=c(min(result[,2]),max(result[,2])),
       ylim=c(min(result[,5]),max(result[,5])))
  m <- apply(cbind(result[,2],result[,5]),2,mean)
  textxy(means,rs,1:a,m=m,cex=1,col="blue")
  origin(m)
  
  
  return(geral.list)
}
