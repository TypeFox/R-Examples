###############################################################################
# Purpose:  Functions for aggregating data from s objects into single         #
#                 dataframes per feature                                      #
# Author:   Ryan Dhindsa                                                      #
###############################################################################

# These functions take data from S object to make dataframes
.write.spike.summary <- function(s){
  # Creates a list of dataframes for spike features. Each df corresponds to 
  #     a single DIV and contains values for each feature
  #
  # Args:
  #   s object
  #           
  # Returns:
  #   list of data frames containing spike data
  
  #initiate empty list to store dataframes
  divs.df = list()
  
  # loop through DIVs in s object and create 1 df per DIV. Each df gets stored in divs.df
  for (i in 1:length(s)) {
    div <- paste("DIV", .get.div(s[[i]]), sep="")
    
    df <- .spike.summary.by.well(s[[i]])
    df = cbind(rownames(df),df)
    df = as.data.frame(unclass(df))  # convert strings to factors
    
    colnames(df)[1] <- "well" 
    divs.df[[div]] = df
  }
  
  divs.df <- do.call(rbind, lapply(names(divs.df), function(x) cbind(div = x, divs.df[[x]])))
  return(divs.df)
}

.compile.ns <-function(s,nspikes) {
  # Calculates network spikes
  # Called from .write.network.spike.summary
  
  active.wells <- .active.wells.network.spikes(nspikes)$ns.all
  if (length(active.wells) >0) {
    # sahar 10292014 - add genotype column and change newcol from 2 to 3
    newcol <-3
    #2 to peak.min and peak.max
    p <- length(active.wells[[1]]$brief) + length(active.wells[[1]]$mean) + newcol
    nsdata <- matrix(0,length(s$well),p)
    temp<-c() # Diana 10/2014
    # Diana Hall 10-31-2014 change
    length.temp.mean<-length(active.wells[[1]]$mean)
    for (j in 1:length(s$well)) {
      cur.well<-s$well[j]
      if ( is.element(cur.well, names(active.wells) ) ){
        temp <- active.wells[[cur.well]]
        nsdata[j,1:length(temp$brief)] <- temp$brief
        nsdata[j,length(temp$brief)+1] <- min(temp$measures[,"peak.val"])
        nsdata[j,length(temp$brief)+2] <- max(temp$measures[,"peak.val"])
        nsdata[j,length(temp$brief)+3] <- s$treatment[cur.well]
        nsdata[j,(length(temp$brief)+newcol+1):p] <- as.double(temp$mean)
        
      } else{
        temp$brief<-c(0, rep(NA,4),0, NA, NA)
        nsdata[j,1:length(temp$brief)] <- temp$brief
        nsdata[j,length(temp$brief)+1] <- NA
        nsdata[j,length(temp$brief)+2] <- NA
        nsdata[j,length(temp$brief)+3] <- NA
        nsdata[j,(length(temp$brief)+newcol+1):p] <- rep(0, length.temp.mean)
      }
    }
    
    nsdata <- data.frame(nsdata)	
    names(nsdata)[1:length(temp$brief)] <- names(active.wells[[1]]$brief)
    names(nsdata)[(length(temp$brief)+1):(length(temp$brief)+newcol)] <- c("peak.min","peak.max","treatment")
    
    for (j in 1:(p-length(temp$brief)-newcol)) {
      names(nsdata)[j+newcol+length(temp$brief)] = paste("t",j,sep = '')
    }
    nsdata<- cbind(s$well,nsdata)
    names(nsdata)[1] <- "well"
    
    return(nsdata)
  }
}
.write.network.spike.summary <- function(s,parameters){
  # Creates a list of dataframes for network spike features. Each df corresponds to 
  #     a single DIV and contains values for each feature
  #
  # Args:
  #   s object
  #           
  # Returns:
  #   list of data frames containing spike data
  
  #initiate empty list to store dataframes
  divs.df = list()
  
  # calculate network spikes
  for (i in 1:length(s)) {
    div <- paste("DIV", .get.div(s[[i]]), sep="")
    
    nspikes.old <- calculate.network.spikes(s[[i]],parameters$sur, parameters$ns.N, parameters$ns.T)
    nspikes <- IGM.summary.network.spikes(s[[i]],nspikes.old,ns.E = 1, parameters$sur)
    basename <- strsplit(basename(s[[i]]$file), "[.]")[[1]][1]
    
    df = .compile.ns(s[[i]],nspikes)
    df = as.data.frame(unclass(df))  # convert strings to factors
    
    df <- df[, -grep("^t[0-9]", colnames(df))]
    
    divs.df[[div]] = df
  }
  
  divs.df <- do.call(rbind, lapply(names(divs.df), function(x) cbind(div = x, divs.df[[x]])))

  return(divs.df)
}

.write.burst.summary <- function(s){
  # Creates a list of dataframes for bursting  features. Each df corresponds to 
  #     a single DIV and contains values for each feature
  #
  # Args:
  #   s object
  #           
  # Returns:
  #   list of data frames containing spike data
  
  masterSum<-.get.burst.info.averaged.over.well(s)
  
  divs.df = list()
  for (i in 1:length(s)){
    div <- paste("DIV", .get.div(s[[i]]), sep="")
    basename <- get.file.basename(s[[i]]$file)
    
    ##########data frame summarized over well
    #get number of object in masterSum[[1]] list
    tempdf<-c(); tempcolnames<-c()
    for (j in 2:length(masterSum[[i]])){
      tempc<-unlist(masterSum[[i]][j])
      tempdf<-cbind(tempdf,tempc)
      tempcolnames<-c(tempcolnames,names(masterSum[[i]][j]) )
    }#end of loop through masterSum list objects
    
    #need to switch around columns so first columns come first
    if (dim(tempdf)[2] > 20) { #for now
      if (dim(tempdf)[1] == 1) {
        df<-cbind(t(tempdf[,21:25]),t(tempdf[,1:20]))
      } else {
        df<-cbind(tempdf[,21:25],tempdf[,1:20])
      }
      colnames<-c(tempcolnames[21:25],tempcolnames[1:20])
      colnames(df)<-colnames
    }
    
    df = as.data.frame(unclass(df))  # convert strings to factors
    df$file = NULL
    row.names(df) = df$well
    divs.df[[div]] = df
  }
  
  divs.df <- do.call(rbind, lapply(names(divs.df), function(x) cbind(div = x, divs.df[[x]])))
  return(divs.df)
}

.create.feat.df <- function(s, df, feature, all_feat_list){
#   feat_df <- dcast(df, well+treatment~treatment, div, value.var = feature)
#   all_feat_list[[feature]] = feat_df
  
  #test
  df=df[!is.na(df$treatment)&df$treatment!="",]
  x = data.frame(df$div, df$well, df$treatment, df[,feature])
  colnames(x) = c("div", "well", "treatment", feature)
  y <- dcast(x, well~div, value.var=feature)
  tempy=cbind(treatment=df[!duplicated(df$well),"treatment"],y[c(-1)])
  y<-cbind(well=y$well,tempy)
  
#   lastdiv = ncol(y) - 1
#   y = y[,c(1,ncol(y), 2:lastdiv)]  # move treatment column from last to second column
   y <- .sort.df(y)
  return(y)
}

.sort.df <- function(df){
  #natural order sorting
  df.divs <- df[3:ncol(df)]
  df.sorted <- df.divs[,mixedorder(names(df.divs)), drop=FALSE]
  df.sorted = cbind(treatment = df$treatment, df.sorted)
  df.sorted = cbind(well = df$well, df.sorted)
  return(df.sorted)
}

IGM.aggregate.feature.data <- function(s, feat.type,parameters=list()){

  # Takes in s object and creates a dataframe for each feature. 
  #     based on the feature type (spikes, ns, etc), it calls appropriate function
  #     (e.g. .write.spike.summary if feat.type is "spikes")
  #
  # Args:
  #   s object
  #   feat.type = "spike", "ns", or "burst"
  #           
  # Returns:
  #   list of data frames (one df per feature)
  
  platename <- get.project.plate.name(s[[1]]$file)
  
  # write feature summaries (calls xxx.summary.by.well from IGM.MEA)
  if (feat.type == "spike") {
    divs.df = .write.spike.summary(s)
  } else if (feat.type == "ns"){
    divs.df = .write.network.spike.summary(s,parameters)
  } else if (feat.type =="burst"){
    divs.df = .write.burst.summary(s)
  }
  

  # create list of dataframes (one dataframe per feature)
  feature.names <- colnames(divs.df)
  remove <- c("div", "treatment", "well")
  feature.names <- setdiff(feature.names, remove)
  
  all.features = list()
  
  #test
  for (i in 1:length(feature.names)) {
    df = .create.feat.df(s, divs.df, feature.names[i], all.features)
    
    all.features[[feature.names[i]]] = df}
#   all.features = sapply(feature.names, function(x) .create.feat.df(divs.df, x, all.features), 
#                         USE.NAMES = FALSE)
  
  
  #all.features = lapply(all.features, function(x) .sort.df(x))
  return(all.features)
}

filter.wells <- function(unfiltered.df, nae,min.electrodes=4){
  # Filters out wells in which there are fewer than 4 active electrodes 
  #    at least 70% of the time
  unfiltered.df=unfiltered.df[!(is.na(unfiltered.df$treatment) | unfiltered.df$treatment==""), ]  # remove wells w/o trt
  
  nae$treatment = NULL
  nae[-1] <- sapply(nae[-1], as.numeric)
  
  num.div <- ncol(nae)-1
  
  inactive <- data.frame(num.inactive = rowSums(nae[, -1] < min.electrodes), total.div = num.div)
  inactive[is.na(inactive$num.inactive),"num.inactive"]=0
  inactive$fraction <- inactive$num.inactive/inactive$total.div
  inactive$well <- nae$well
  
  active.wells <- with(inactive, {subset.data.frame(inactive, fraction < 0.7, select = well)})
  
  filtered.df = unfiltered.df[unfiltered.df$well %in% active.wells$well, ]
  
  #replace na's with 0's
  filtered.df[filtered.df == "NaN"] = NA  # first replace NaN with NA
  filtered.df[is.na(filtered.df)] <- 0    # then replace NA's with 0
  
  return(filtered.df)
}

write.features.to.files <- function(s, features.list, output.dir, type) {
  # Takes in list of dataframes (one per feature) and writes out each 
  #     df to a csv file
  #
  # Args:
  #   s object
  #   features.list = list of dataframes
  #   output.dir = directory where files will be put (will make separate folders
  #                  ns, spikes, and bursts) 
  #           
  # Returns:
  #   one csv per feature
  
  # change to create subdir for each file type
  platename <- get.project.plate.name(s[[1]]$file)
  out.folder <- paste0(output.dir,"/",type)
  dir.create(out.folder, showWarnings = FALSE)
  invisible(sapply(names(features.list), 
                   function (x) write.csv(features.list[[x]], file=paste0(out.folder, "/", platename, "_", x, ".csv"), 
                                          row.names = F)))
}

