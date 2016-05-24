# target and antibody should already be combined to 'target', slide and sixwell combined to ID
getBlocks <- function(subsample, time=0, ab=NULL, distinguish) {
  
  mysub <- subset(subsample, time==time & ab == ab)
  N <- dim(mysub)[1]
  IDs <- unique(mysub$ID)
  #stimulations <- unique(paste(mysub$treatment, mysub$cellline))
  stimulations <- unique(do.call(paste, mysub[,distinguish]))
  
  M <- matrix(0, ncol=length(IDs), nrow=length(stimulations))
  
  for(i in 1:N) {
    #myrow <- which(stimulations == paste(mysub$treatment[i], mysub$cellline[i]))
    myrow <- which(stimulations == do.call(paste, mysub[,distinguish])[i])
    mycol <- which(IDs == mysub$ID[i])
    M[myrow, mycol] <- 1
  }
  
  
  ## Analyze matrix with zero and one entries for connection components
  analyzeBlocks <- function(M) {
    out <- which(apply(M, 1, sum)==0)
    if(length(out) > 0) {
      M <- M[-out,]
      cat("matrix contains zero rows which have been eliminated\n")
    }
    n <- dim(M)[1]
    rcomponents <- list()
    ccomponents <- list()
    
    counter <- 0
    while(length(unlist(rcomponents))< n) {
      
      counter <- counter + 1
      
      if(length(unlist(rcomponents)) == 0) {
        w <- 1
      } else {
        mysample <- (1:n)[-unlist(rcomponents)]
        w <- min(mysample)
      }
      
      repeat {	
        
        v <- unique(rapply(as.list(w), function(i) which(M[i,] == 1)))
        wnew <- unique(rapply(as.list(v), function(j) which(M[,j] == 1)))
        if(length(wnew) == length(w)) break
        w <- wnew
      }
      rcomponents[[counter]] <- w
      ccomponents[[counter]] <- v
      
      
    }
    
    #print(M[unlist(rcomponents), unlist(ccomponents)])
    
    
    return(rcomponents)   
  
  }
   
  return(analyzeBlocks(M))
  
}



averageData <- function(subsample, scaling = c("slide", "replicate"), distinguish = c("cellline", "treatment")) {
  
  
  IDs <- unique(do.call(paste, subsample[,scaling]))
  print(IDs)
  subsample$ID <- match(do.call(paste, subsample[,scaling]), IDs)
  
  result <- c()
  targets <- unique(subsample$ab)
  for(A0 in targets) {
    mytarget <- which(targets==A0)
    cat(paste("starting nlsFit for target ", A0, "(",mytarget,"/",length(targets),")\n",sep=""))
    ab<-NULL
    subsample0 <- subset(subsample, ab == A0)
    cat("analyzing independent stimulations ...")
    indepStim <- getBlocks(subsample0, ab=A0, distinguish=distinguish)
    cat(paste(" done. Found", length(indepStim),"subsets\n"))
    
    stimulations <- unique(do.call(paste, subsample0[,distinguish]))
    
    
    for(A2 in 1:length(indepStim)) {
      #subsample1 <- subset(subsample0, paste(treatment, cellline) %in% stimulations[indepStim[[A2]]])
      subsample1 <- subsample0[do.call(paste, subsample0[,distinguish]) %in% stimulations[indepStim[[A2]]], ]
      
      nTimes <- length(unique(subsample1$time))
      newsignals <- c()
      cat("producing new signal data frame ...")
      for(A1 in unique(subsample1$time)) {
        subsample2 <- subset(subsample1, time==A1)
        
        
        #IDs <- unique(paste(subsample2$ID, subsample2$treatment, subsample2$cellline))
        IDs <- unique(do.call(paste, subsample2[,c("ID", distinguish)]))
        #print(rev(IDs)[1])
        #print(IDs)
        signals <- c()
        for(IDvalue in IDs) {
          
          #alldata <- subset(subsample2, paste(ID, treatment, cellline)==IDvalue)
          alldata <- subsample2[do.call(paste, subsample2[,c("ID", distinguish)])==IDvalue, ]
          #print(alldata)
          #print(dim(alldata))
          meanSignal <- mean(alldata$signal)
          varSignal <- var(alldata$signal)
          var0 <- mean(alldata$var0)
          varR <- mean(alldata$varR)
          myID <- unique(alldata$ID)
          mystimulation <- unique(alldata$treatment)
          myzelllinie  <- unique(alldata$cellline)
          
          sig1 <- data.frame(meanSignal = meanSignal, varSignal = varSignal, ID = myID, var0 = var0, varR=varR)
          sig2 <- cbind(sig1, alldata[1,distinguish])
          colnames(sig2) <- c(colnames(sig1), distinguish)
          
          signals <- rbind(signals, sig2)
          
        }
        
        newsignals <- rbind(newsignals, cbind(data.frame(time=A1), signals))
              
      }
      #newsignals <- newsignals[order(newsignals$time, newsignals$treatment, newsignals$cellline, newsignals$ID),]
      newsignals <- newsignals[do.call(order, newsignals[,c("time", distinguish, "ID")]),]
      cat(" done\n")
      
      
      response <- newsignals$meanSignal
      zeros <- rep(0, length(response))
      variances <- newsignals$varSignal
      var0 <- newsignals$var0
      varR <- newsignals$varR
      
      #ys <- unique(paste(newsignals$time, newsignals$treatment, newsignals$cellline))
      ys <- unique(do.call(paste, newsignals[,c("time", distinguish)]))
      #indlabels <- sapply(as.list(1:length(ys)), function(i) which(paste(newsignals$time, newsignals$treatment, newsignals$cellline) == ys[i])[1])
      indlabels <- sapply(as.list(1:length(ys)), function(i) which(do.call(paste, newsignals[, c("time", distinguish)]) == ys[i])[1])
      ss <- unique(newsignals$ID)
      N <- dim(newsignals)[1]
      
      indy <- c()
      inds <- c()
      
      for(i in 1:N) {
        #indy <- c(indy, which(ys == paste(newsignals[i,"time"], newsignals[i,"treatment"], newsignals[i,"cellline"])))
        indy <- c(indy, which(ys == do.call(paste, newsignals[,c("time", distinguish)])[i]))
        inds <- c(inds, which(ss == newsignals[i, "ID"]))
      }
      
      model <- function(y, s0) {
        s <- c(1, s0)
        
        
        yisj <- sapply(as.list(1:N), function(i) y[indy[i]]/s[inds[i]])
        
        
        mij <- function(y,s) {
          mij <- (yisj - response)/(sqrt(var0+varR*yisj^2))
          return(mij)
        }
        
        
        grad <- t(sapply(as.list(1:N), function(i) {
          resy <- rep(0, length(y))
          ress <- rep(0, length(s))
          resy[indy[i]] <- ((var0[i]+varR[i]*yisj[i]*response[i])/((var0[i] + varR[i]*yisj[i]^2)^(1.5)))*(1/s[inds[i]])
          ress[inds[i]] <- -((var0[i]+varR[i]*yisj[i]*response[i])/((var0[i] + varR[i]*yisj[i]^2)^(1.5)))*yisj[i]/s[inds[i]]
          return(c(resy, ress[-1]))
        }))
        
          
        m <- mij(y,s)
        
        attr(m, "gradient") <- grad 
        return(m)
        
      }
      #ystart <- sapply(as.list(unique(paste(newsignals$time, newsignals$treatment, newsignals$cellline))), function(z) mean(subset(newsignals, paste(time, treatment, cellline)==z)$meanSignal))
      ystart <- sapply(as.list(unique(do.call(paste, newsignals[,c("time", distinguish)]))), function(z) mean(newsignals[do.call(paste, newsignals[, c("time", distinguish)])==z, "meanSignal"]))
      
      sstart <- rep(1, -1 + length(unique(newsignals$ID)))
      
      cat("fitting true signals and scaling factors ...")
      
      fit <- try(nls(zeros ~ model(y,s), start=list(y = ystart, s=sstart), control=list(warnOnly=TRUE)), silent=FALSE)
      
      if(class(fit)=="try-error") {
        cat(" unable to fit parameters\n")
      } else {
        cat(" done\n")
        
        y <- coef(fit)[1:length(ystart)]
        sy <-sqrt(diag(vcov(fit)))[1:length(ystart)] 
        
        res1 <- data.frame(ID = newsignals$ID[indlabels], connection=A2, ab=A0, time=newsignals$time[indlabels], signal=y/mean(y), sigma=sy/mean(y))
        res2 <- cbind(newsignals[indlabels, distinguish], res1)
        colnames(res2) <- c(distinguish, colnames(res1))
        
        
        #result <- rbind(result, data.frame(cellline = newsignals$cellline[indlabels], treatment = newsignals$treatment[indlabels], ID = newsignals$ID[indlabels], connection=A2, ab=A0, time=newsignals$time[indlabels], signal=y/mean(y), sigma=sy/mean(y)))
        result <- rbind(result, res2)
      }
      
      
      
    }
    
  }
  
  return(result)
  
  
}
