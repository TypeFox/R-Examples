thsu <-
function(data,interaction=FALSE){  
    a <- nrow(data)
    b <- ncol(data)
    data.m <- as.matrix(data)
    l.data <- length (data.m)
    data.r <- matrix(NA,a,b) 
    ranks <- matrix(NA,a,b)
    ranks.y <- matrix(NA,a,b)
    np.1 <- matrix(NA,a,b)
    np.2 <- matrix(NA,a,b)
    np.3 <- matrix(NA,a,b)
    np.4 <- matrix(NA,a,b)    
    np1 <- numeric()
    np2 <- numeric()
    np3 <- numeric()
    np4 <- numeric()
    k <- 2
    
    for (i in 1:nrow(data.m)){
      for (j in 1:ncol(data.m)){
        data.r[i,j] <- (data.m[i,j]) - (mean(data.m[i,]))
      }
    }
    
    for (j in 1:ncol(data.r)){
      
      ranks[,j] <- rank(-data.r[,j])
    }
    
    for (j in 1:ncol(data.m)){
      
      ranks.y[,j] <- rank(-data.m[,j])
    }
    
    for (i in 1:nrow(data)){
      for (j in 1:ncol(data)){
        
        np.1[i,j] <- abs((ranks[i,j])-median(ranks[i,]))
        np.2[i,j] <- abs((ranks[i,j])-median(ranks[i,]))/(median(ranks.y[i,]))
        np.3[i,j] <- ((ranks[i,j]-mean(ranks[i,]))^2) / b
      }
      np1[i] <- round((1/b) * (sum(np.1[i,])),digits=4)
      np2[i] <- round((1/b) * (sum(np.2[i,])),digits=4)
      np3[i] <- round((sqrt(sum(np.3[i,]))) / (mean(ranks.y[i,])),digits=4)
    }
    
    for (i in 1:nrow(data)){
      for (j in 1:(b-1)){
        
        np.4[i,j] <- abs(ranks[i,j] - ranks[i,k]) 
        while(k < b)
          k <- k + 1
      }
      np4[i] <- round((2/(b*(b-1))) * (sum((np.4[i,j]) / (mean(ranks.y[i,])))),digits=4) 
    }
    
    means <- round(as.numeric(rowMeans(data)),digits=4)
    result <- as.data.frame(cbind(rownames(data),means,np1,np2,np3,np4))
    colnames(result) <- c("Gen","Mean","N1","N2","N3","N4")
    
    rank.y <- apply(-data,2,rank)
    ranks.sum.y <- apply(rank.y,1,sum)
    sd.rank = round(apply(rank.y,1,sd),digits=4)
    ranks.y = data.frame(rank.y,ranks.sum.y,sd.rank)
    colnames(ranks.y) = c(colnames(rank.y),"Sum", "Sd")
    cor.rank.y <- round(cor(ranks.y, method="pearson"), digits = 4)
    geral.list <- list("ThSu"=result,"Ranks"=ranks.y,"Correlations"=cor.rank.y)
    
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
    
    n1 <- as.numeric(result[,3])
    n2 <- as.numeric(result[,4])
    n3 <- as.numeric(result[,5])
    n4 <- as.numeric(result[,6])
    
    plot(means,n1,pch=19,cex=0.5,main="Means x N1",xlab=expression(Mean[Phenotypic]),
         ylab=expression(N[1]),xlim=c(min(means),max(means)),
         ylim=c(min(n1),max(n1)))
    m <- apply(cbind(means,n1),2,mean)
    textxy(means,n1,1:a,m=m,cex=1,col="blue")
    origin(m)
    plot(means,n2,pch=19,cex=0.5,main="Means x N2",xlab=expression(Mean[Phenotypic]),
         ylab=expression(N[2]),xlim=c(min(means),max(means)),
         ylim=c(min(n2),max(n2)))
    m <- apply(cbind(means,n2),2,mean)
    textxy(means,n2,1:a,m=m,cex=1,col="blue")
    origin(m)
    plot(means,n3,pch=19,cex=0.5,main="Means x N3",xlab=expression(Mean[Phenotypic]),
         ylab=expression(N[3]),xlim=c(min(means),max(means)),
         ylim=c(min(n3),max(n3)))
    m <- apply(cbind(means,n3),2,mean)
    textxy(means,n3,1:a,m=m,cex=1,col="blue")
    origin(m)
    plot(means,n4,pch=19,cex=0.5,main="Means x N4",xlab=expression(Mean[Phenotypic]),
         ylab=expression(N[4]),xlim=c(min(means),max(means)),
         ylim=c(min(n4),max(n4)))
    m <- apply(cbind(means,n4),2,mean)
    textxy(means,n4,1:a,m=m,cex=1,col="blue")
    origin(m)
    
    
    return(geral.list)
  }
