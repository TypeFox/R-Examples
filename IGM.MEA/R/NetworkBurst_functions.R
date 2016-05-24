# convert the spike list to a matrix, probably not necessary anymore
# but just to be consistent with matlab for now
.NB.prepare <- function(S) {
  spikes <- S$spikes
  colnames <- names(spikes)
  min.values <- sapply(spikes,min)
  max.values <- sapply(spikes,max)
  
  slength = sapply(spikes,length)
  max.elements <- max(slength)
  
  data <- matrix(0,length(colnames),max.elements)
  for (i in  1:length(colnames)) {
    data[i,1:slength[[i]]] <- spikes[[i]]
  }
  rownames(data) <- colnames
  data<- t(data)
  data
}

#generate binned data for the well
.NB.select.and.bin <- function(data, well,sbegin,send,bin) {
  well.data <- data[,grep(well,colnames(data))]
  well.data[well.data<sbegin | well.data>send] <- 0
  if (is.null(dim(well.data))) {
    dim(well.data) <- c(length(well.data),1) 
  }
  well.data <- well.data[rowSums(well.data)>0,]
  well.data <- round(well.data * 12500)
  temp <- well.data[well.data>0]
  if (length(temp)>0) {
    t_min <- min(temp)
    t_max <- max(temp)
    well.data <- well.data - t_min + 1
    range <- (t_max-t_min+1)
    bins <- floor(range / bin)
    if (range %% bin > 0) {
      bins <- bins + 1
    }
    x <- (0:bins)*bin + 1
    y <- c(-Inf,x, max(x)+bin, +Inf)
    
    if (is.null(dim(well.data))) {
      dim(well.data) <- c(length(well.data),1) 
    }
    well.data.new <- matrix(0,bins+1,dim(well.data)[2])
    for (i in 1: dim(well.data)[2]) {
      temp <- hist(well.data[,i],breaks = y, plot = FALSE)$counts
      temp <- temp[2:(length(temp)-1)]
      well.data.new[,i] <- temp
    }
    well.data.new[well.data.new>1] <-1
  } else {
    well.data.new <- temp
  }
  well.data.new
}

.NB.Gaussian.Filter <- function(sigma, well.data, min_electrodes) {
  sigma.input <- sigma
  if (length(well.data) == 0) {
    f2 <- numeric()
  } else {
    #a very flat filter. Consider a much shapper one later
    if (sigma %% 2 == 0) { sigma <- sigma+1}
    filter_size <- sigma
    half_size <- (filter_size-1) / 2
    x <- seq(-half_size,half_size,length=filter_size)
    gaussFilter <- exp(-x ^ 2 / (2 * sigma ^ 2))
    gaussFilter <- gaussFilter / sum (gaussFilter)
    f <- well.data
    if (length(which(apply(well.data, 2, max)>0)) >= min_electrodes) {
      for (i in 1:dim(f)[2]) {
        if (max(well.data[,i]) > 0) {
          f[,i] <- filter(well.data[,i],gaussFilter,circular = TRUE)
          #temp <- filter(well.data[,i],gaussFilter,circular = TRUE)
          #f[,i] <- temp[(half_size+1):(length(temp) - half_size)]
          
          f[,i] <- f[,i] / max(f[,i])
        }
      }
      f1 <- rowSums(f)
      f1 <- f1/max(f1)
      f2 <- filter(f1,gaussFilter,circular = TRUE)
      #f2 <-conv(gaussFilter,f1)
      #f2 <- f2[(half_size+1):(length(f2) - half_size)]
      f2 <- f2/ max(f2)
      dim(f2) <- c(length(f2),1)
      colnames(f2) <- sigma.input
    } else {
      f2 <- numeric()
    }
  }
  f2
}

.NB.Get.Spike.Stats <- function(well.data,timespan) {
  stat <- matrix(0,1,3)
  if (length(well.data) > 0) {
    stat[1] = length(which(apply(well.data, 2, max)>0)) #number of active electrodes
    stat[2] = sum(well.data)/timespan; #well level mfr
    if (stat[1] > 0)  {     #mfr by active electrodes
      stat[3] <- stat[2] / stat[1]
    }
  }
  colnames(stat) <- c("n.active.electrodes","mfr.per.well","mfr.per.active.electode")
  stat
}

.NB.Get.Burst.Stats.Intensities <- function(well.data,f,timespan,bin.time,min_electrodes) {
  n <- length(f)
  stat <- matrix(-1,11,n)
  rownames(stat) <- c("mean.network.burst.time.per.second",
                      "number.of.spikes.in.network.bursts",
                      "percentage.of.spikes.in.network.bursts",
                      "spike.intensity",
                      "spike.intensity.by.electrodes",
                      "number.of.spikes.in.network.bursts",
                      "number.of.spikes.per.network.burst",
                      "mean.spikes.per.active.electrode.in.network.burest",
                      "mean.spikes.per.active.electrode.in.network.burest.per.second",
                      "mean.spikes.per.well.in.network.burest",
                      "total.number.of.network.bursts")
  colnames(stat) <- rep("NA",n)
  for (current in 1:n) {
    if (length(f[[current]] > 0)) { #not all have names
      colnames(stat)[current] <- colnames(f[[current]]) # not all have names 
    } 
  }
  
  for (current  in 1:n) {
    temp = f[[current]]
    if (length(temp)>0) { 
      nActiveElectrodes <- length(which(apply(well.data, 2, max)>0))
      temp[temp<0] <- 0
      level <- .graythresh(temp)
      indicator0 <- temp >= level #get timestamps that are identified as part of network bursts
      
      if (min_electrodes > 0) {
        indicator = c(0,indicator0,0) #prepare for finding intervals
        diff_indicator <- diff(indicator)
        burst_start <- which(diff_indicator==1)
        burst_end <- which(diff_indicator==-1)
        if (length(burst_start) != length(burst_end)) { 
          stop('Burst region no match')
        }
        
        ##filtered regions based on minimum number of active electrodes
        for (region in 1: length(burst_start)) {
          current_region_nAE <- length(which(colSums(well.data[(burst_start[region]):(burst_end[region]-1),])>0))
          if (current_region_nAE < min_electrodes) {
            indicator0[burst_start[region]:(burst_end[region]-1)] <- 0
          }
        }
      }
      
      indicator = c(0,indicator0,0) #prepare for finding intervals
      diff_indicator <- diff(indicator)
      burst_start <- which(diff_indicator==1)
      burst_end <- which(diff_indicator==-1)
      if (length(burst_start) != length(burst_end)) {
        stop('Burst region no match')
      }
      
      stat[1,current] <- sum(indicator0) * bin.time / timespan #mean network burst time per second
      totalspikes <- rowSums(well.data)
      SpikesInNetworkBursts <- sum(totalspikes*indicator0)
      stat[2,current] <- SpikesInNetworkBursts
      stat[3,current] <- stat[2,current] / sum(totalspikes) #percentage of spikes in network bursts
      
      TotalBurstRegions <- sum(burst_end - burst_start)
      stat[4,current] <- SpikesInNetworkBursts / TotalBurstRegions #overall Spike Intensity 
      stat[5,current] <- stat[4,current] / nActiveElectrodes
      stat[6,current] <- SpikesInNetworkBursts / length(burst_start)
      stat[7,current] = stat[6,current] / nActiveElectrodes
      
      burst_region_length <- burst_end - burst_start
      total_spikes_per_active_electrode <- totalspikes / nActiveElectrodes
      spikes_in_region_per_active_electrode <- matrix(0,length(burst_start),2)
      for (region in 1: length(burst_start)) {
        spikes_in_region_per_active_electrode[region,1] <- sum(total_spikes_per_active_electrode[burst_start[region]:(burst_end[region]-1)])
        #normalized to HZ per active electrode within burst region
        spikes_in_region_per_active_electrode[region,2] <- spikes_in_region_per_active_electrode[region,1] / (burst_region_length[region] * bin.time)
      }
      stat[8,current] <- mean(spikes_in_region_per_active_electrode[,1])
      stat[9,current] <- mean(spikes_in_region_per_active_electrode[,2])
      
      stat[10,current] = stat[8,current]*nActiveElectrodes
      stat[11,current] = length(burst_start)
    }
  }  
  stat
}

#do not change bin and sf for MEA recordings, skip is in seconds
.NB.Extract.Fetures <- function(S, Sigma, min_electrodes, local_region_min_nAE,duration = 0, bin = 25,sf = 12500,skip = 10) {
  df.spikes <- .NB.prepare(S)
  wells <- unique(substr(colnames(df.spikes),1,2))
  output <- list()
  #time interval
  #skip <- 10 #10 seconds
  sbegin <- floor(min(df.spikes[df.spikes>0]/skip))* skip + skip
  send <- floor(max(df.spikes[df.spikes>0]/skip))* skip
  timespan <- send - sbegin
  #bin = 25; #2ms, this is a parameter based on our recording
  bin.time <- bin / sf
  
  for (well.index in 1: length(wells)) {
    #print(well.index)  
    well <- wells[well.index]
    well.data <- .NB.select.and.bin(df.spikes, well,sbegin,send,bin)	 
    f <- list()
    for (i  in 1: length(Sigma)) {
      f[[i]] <- .NB.Gaussian.Filter(Sigma[i], well.data, min_electrodes)
    }
    stat0 <- .NB.Get.Spike.Stats(well.data,timespan)	
    stat1 <-.NB.Get.Burst.Stats.Intensities(well.data,f,timespan,bin.time,local_region_min_nAE)
    
    output[[well.index]] <- list()
    output[[well.index]]$stat0 <- stat0
    output[[well.index]]$stat1 <- stat1
    output[[well.index]]$well <- well
    
  }
  result <- list(data = output, DIV = unlist(strsplit(unlist(strsplit(basename(S$file),"_"))[4],"[.]"))[1])
  result
}

.NB.Merge.Result <- function(s,result,Sigma) {
  Wells <- character()
  DIVs <- character()
  Phenotypes <- character()
  Data <- numeric()
  
  Count <- 0
  w.fun <- function(x) {x$well}
  for (i in 1: length(result)) {
    current <- result[[i]]  
    current.Wells <- sapply(current$data, w.fun)
    DIVs <- c(DIVs, rep(current$DIV,length(current.Wells)))
    Wells <- c(Wells,current.Wells)
    Phenotypes<- c(Phenotypes, s[[i]]$treatment[current.Wells])
    
    data <- matrix(0,length(current.Wells),length(c(as.vector(current$data[[1]]$stat1),as.vector(current$data[[1]]$stat0))))
    for (j in 1:length(current.Wells)) {
      data[j,] <- c(as.vector(current$data[[j]]$stat1),as.vector(current$data[[j]]$stat0))
    }
    Data <- rbind(Data,data)
  }
  
  
  feature.names <- character()
  for (i in 1:length(Sigma)) {
    feature.names <- c(feature.names,paste(rownames(current$data[[1]]$stat1), Sigma[i], sep = "_"))
  }
  feature.names <- c(feature.names, colnames(current$data[[1]]$stat0))
  
  df = data.frame(DIVs, Wells, Phenotypes,Data)
  names(df)[4:dim(df)[2]] <- feature.names
  df
}

NB.matrix.to.feature.dfs <- function(data) {
  n.features = dim(data)[2] - 3 #escape the first columns
  dfs <- list()
  data[,'Wells'] <- factor(data[,'Wells']) #drop unused levels
  Well.stat <- table(data[,'Wells'])
  
  ref.matrix <- matrix(NA,length(Well.stat),max(Well.stat))
  Wells <- unique(data[,c('Wells','Phenotypes')])
  rownames(ref.matrix) <- Wells[,'Wells']
  colnames(ref.matrix) <- unique(data[,1])
  
  #now change columnnames to match Ryan's code
  colnames(Wells) <- c('well','treatment')
  n <- dim(data)[1]
  for (index in 1:n.features) {
    data.matrix <- ref.matrix
    for (i in 1:n) {
      data.matrix[data[i,'Wells'],data[i,'DIVs']] <- data[i,index + 3]
    }
    dfs[[index]] <- cbind(Wells,data.matrix)
  }
  names(dfs) <- colnames(data)[4:(n.features+3)]
  dfs
}

.wilcox.test.perm <- function(data,np,g1, g2, feature.index) {
  #now figure out the p from data
  d1 <- data[data[,'Phenotypes']==g1,feature.index]
  d1 <- d1[d1>=0]
  d2 <- data[data[,'Phenotypes']==g2,feature.index]
  d2 <- d2[d2>=0]
  suppressWarnings(data.p <- wilcox.test(d1, d2)$p.value)
  
  #subsetting data to genotypes and feature, and also reformat into matrix for easy permutation
  d <- data[(data[,'Phenotypes']==g1) | (data[,'Phenotypes']==g2),c(2,feature.index)]
  d[,'Wells'] <- factor(d[,'Wells']) #drop unused levels
  Well.stat <- table(d[,'Wells'])
  data.matrix <- matrix(-1,length(Well.stat),max(Well.stat))
  rownames(data.matrix) <- names(Well.stat)
  n.cases <- length(unique(data[data[,'Phenotypes']==g1,'Wells']))
  n <- dim(data.matrix)[1]
  for (i in 1:n) {
    temp <- d[d[,'Wells'] == rownames(data.matrix)[i],2]
    data.matrix[i,1:length(temp)] <- temp
  }
  
  outp <- matrix(0,np,1)
  for (i in 1:np) {
    cases = sample(n,n.cases)
    d1 = as.vector(data.matrix[cases,])
    d1 <- d1[d1>=0]
    d2 = as.vector(data.matrix[-cases,])
    d2 <- d2[d2>=0]
    suppressWarnings(outp[i] <- wilcox.test(d1, d2)$p.value)
  }
  outp = sort(outp)
  
  perm.p <- length(which(outp<data.p)) / np
  result <- list(perm.p = perm.p, outp = outp)
  result
}

calculate.network.bursts <- function(s,Sigma,min_electrodes,local_region_min_nAE) {
  # extract features and merge features from different recordings into one data frame
  R <- NULL
  if (length(s)>0) {
    result <- list()
    for (i in 1:length(s)) {
      result[[i]] <- .NB.Extract.Fetures(s[[i]], Sigma, min_electrodes, local_region_min_nAE) 
    }
    R <- .NB.Merge.Result(s,result,Sigma)
  }
  R
}
