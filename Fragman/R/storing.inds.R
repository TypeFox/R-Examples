storing.inds <-
  function (folder, channels=NULL, fourier=TRUE, saturated=TRUE, lets.pullup=TRUE, plotting=FALSE, rawPlot=FALSE) {
    #if (!require("seqinr")) {
    #  install.packages("seqinr")
    #  require("seqinr")
    #}
    #storing the names of the files
    coli <-  c("cornflowerblue", "chartreuse4", "gold2", "red", 
               "orange", "purple")
    setwd(paste(folder))
    listp2 <- dir(folder, "*.fsa")
    # list to store all matrices of 4 columns
    all.inds.mats <- list(NA)
    
    print("Reading FSA files")
    
    ####################
    ## initialize the progress bar
    count <- 0
    tot <- length(listp2)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
    #####################
    
    for(i in 1:length(listp2)){
      ###################
      count <- count + 1
      ###################
      fsaFile<-read.abif(listp2[i]) # file read
      lens <- lapply(fsaFile$Data, length)
      aaa <- table(unlist(lens)) # length of vectors in data
      #### CONDITION IF CHANNELS WAS NOT SPECIFIED
      if(is.null(channels)){
        channels <- as.vector(aaa[which(as.numeric(names(aaa)) > 1000)])
      }
      real.len <- as.numeric(names(aaa)[which(aaa == channels & as.numeric(names(aaa)) > 1000 )]) # length of the data elements that meet the requirement to be a data file
      v <- as.vector(which(unlist(lens) == real.len)) # length to find, elements storing the data found
      reads<-list(NA) # to store info
      for(j in 1:channels){ # for each channel
        v2 <- v[j]
        reads[[j]] <- fsaFile$Data[[v2]]
      }
      all.inds.mats[[i]] <- matrix(unlist(reads),ncol=channels)
      names(all.inds.mats)[i] <- as.character(listp2[i])
      ################################
      setTxtProgressBar(pb, count/tot)### keep filling the progress bar
      ################################
    }
    close(pb) # close the progress bar
    if(fourier == TRUE){ #FOURIER
      print("Applying Fourier tranformation for smoothing...")
      all.inds.mats <- lapply_pb(all.inds.mats, function(x){apply(x, 2, transfft)})
    }
    if(saturated == TRUE){ #SATURATION
      print("Checking and correcting for saturated peaks...")
      all.inds.mats <- lapply_pb(all.inds.mats, function(x){apply(x, 2, saturate)})
    }
    if(lets.pullup==TRUE){#PULL UP
      print("Applying pull up correction to the samples to decrease noise from channel to channel")
      if(plotting == TRUE){
        all.inds.mats <- lapply_pb(all.inds.mats, pullup, channel=channels, plotting=TRUE)
      }
      all.inds.mats <- lapply_pb(all.inds.mats, pullup, channel=channels)
    }
    ### ------------------
    ### ------------------
    
    if(rawPlot == TRUE){
      layout(matrix(1:2,2,1))
      coli <- cfp <- c("cornflowerblue", "chartreuse4", "gold2", "red", "orange", "purple")
      naname <- c("FAM","HEX","NED","ROX","LIZ")
      
      print("Plotting raw data")
      
      ####################
      ## initialize the progress bar
      count <- 0
      tot <- length(listp2)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
      #####################
      for(i in 1:dim(all.inds.mats[[1]])[2]){
        plot(all.inds.mats[[1]][,i], col=transp(coli[i],0.6), type="l", ylab="RFU", main=paste(naname[i],"channel. Plot",i,"of",dim(all.inds.mats[[1]])[2]), las=2)
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey75")
        
        if(length(all.inds.mats) > 1){
          for(j in 1:length(all.inds.mats)){
            ###################
            count <- count + (1/dim(all.inds.mats[[1]])[2])
            ###################
            lines(all.inds.mats[[j]][,i], col=transp(coli[i],0.4), lwd=.2)
            ################################
            setTxtProgressBar(pb, count/tot)### keep filling the progress bar
            ################################
          }
        }
      }
      close(pb) # close the progress bar
    }
    
  
    layout(matrix(1,1,1))
    return(all.inds.mats)
  }

