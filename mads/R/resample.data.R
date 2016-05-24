#' Resamples the data for the bootstrap
#'
#' Can either resample samples within strata or observations within samples
#'
#' Once the observations have been resampled duplicates are given a different
#' object id and the obs.table is updated appropriately. 
#'  
#' @param resample if "samples" resample samples within strata, if 
#'   "observations" resample observations within samples.
#' @param obs.table dataframe of observation records with fields object,
#'   Region.Label, and Sample.Label which give links to sample.table,
#'   region.table and the data records used in \code{model}
#' @param sample.table dataframe of sample records - Region.Label,
#'   Sample.Label, Effort
#' @param ddf.dat.master list of complete/original datasets
#' @param double.observer boolean indicating if it is a double observer survey
#' @param subset.variable character describing which variable in the dataset
#'   should be used in subsetting the data.
#' @return a list with 2 elements: 
#'   ddf.dat.working a list of resampled datasets to be used in the analyses 
#'   obs.table an updated obs.table with additional entries for data replicates
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @keywords data manipulation
#'
resample.data <- function(resample, obs.table, sample.table, ddf.dat.master, double.observer, subset.variable = "species"){
# resample.data function to non-parametrically resample the observations
#
# Arguments:
#   resample       - character option either "samples" or "observations"
#   obs.table      - dataframe of observation records
#   sample.table   - dataframe of sample records
#   ddf.dat.master - list of original datasets
#
# Value: list containind the following elements
#   ssf.dat.working - a list of resampled datasets
#   obs.table       - an updated obs.table
#
# Functions Used: renumber.duplicates 
#
  #create data storage
  ddf.dat.working <- list()
  species.name <- names(ddf.dat.master)

  #RESAMPLE SAMPLES WITHIN STRATA
  if(resample == "samples"){
    #Test to make sure that region labels are unique
    #if(length(unique(region.dat$Region.Label))!=length(region.dat$Region.Label))
    #  stop("Region labels must be unique")
    #samples=merge(region.dat,sample.dat,by.x="Region.Label",all.x=TRUE,all.y=TRUE)
    region.names <- as.character(unique(sample.table$Region.Label))
    
    #for every strata
    for(r in seq(along = region.names)){
      #get sample IDs
      temp <- sample.table[sample.table$Region.Label == region.names[r],]
      #Test to see if all sample labels are unique within strata?
      #resample
      new.samples <- data.frame(Sample.Label = sort(sample(temp$Sample.Label, nrow(temp), replace = TRUE)), new.Sample.Label = rep(NA,nrow(temp)))
      new.samples$Sample.Label <- as.character(new.samples$Sample.Label)
      new.samples$new.Sample.Label[1] <- new.samples$Sample.Label[1]
      #Rename duplicate samples
      count <- 1
      for(samp in seq(along = new.samples$Sample.Label)[-1]){
        if(new.samples$Sample.Label[samp] != new.samples$Sample.Label[samp-1]){
          new.samples$new.Sample.Label[samp] <- new.samples$Sample.Label[samp]
          count <- 1  
        }else{
          new.samples$new.Sample.Label[samp] <- paste(new.samples$Sample.Label[samp],"rep",count, sep = "") 
          count <- count + 1
        }
      }
      new.samples$Sample.Label <- as.factor(new.samples$Sample.Label)
      new.samples$new.Sample.Label <- as.factor(new.samples$new.Sample.Label)
      #get observation numbers
      if(r == 1){
        resample.result <- merge(new.samples, obs.table, by = "Sample.Label")
        resample.result <- resample.result[order(resample.result$object),]
        new.sample.table <- merge(new.samples, sample.table, by = "Sample.Label")
      }else{
        resample.result <- rbind(resample.result, merge(new.samples, obs.table, by = "Sample.Label"))
        new.sample.table <- rbind(new.sample.table, merge(new.samples, sample.table, by = "Sample.Label"))
      }
    }
    #get data for chosen observations 
    for(sp in seq(along = species.name)){
      #need to add data sample by sample so that the new sample ID's can be added
      new.samples <- new.sample.table$new.Sample.Label
      temp.data <- merge(ddf.dat.master[[species.name[sp]]], resample.result[resample.result$new.Sample.Label == new.samples[1],], by='object')[,1:ncol(ddf.dat.master[[species.name[sp]]])]
      ddf.dat.working[[species.name[sp]]] <- cbind(temp.data, new.sample.id = rep(new.samples[1], nrow(temp.data)))
      for(samp in seq(along = new.samples)[-1]){
        temp.data <- merge(ddf.dat.master[[species.name[sp]]], resample.result[resample.result$new.Sample.Label == new.samples[samp],], by='object')[,1:ncol(ddf.dat.master[[species.name[sp]]])]
        temp.data <- cbind(temp.data, new.sample.id = rep(new.samples[samp], nrow(temp.data)))
        ddf.dat.working[[species.name[sp]]] <- rbind(ddf.dat.working[[species.name[sp]]], temp.data)
      }
      #renumber duplicates and update obs.table
      renumber.duplicate.results <- renumber.duplicates(ddf.dat.working[[species.name[sp]]], obs.table, double.observer)
      ddf.dat.working[[species.name[sp]]] <- renumber.duplicate.results[[1]]
      obs.table <- renumber.duplicate.results[[2]]
      #Store reference from object id to new sample id
      if(sp == 1){
        all.data <- ddf.dat.working[[species.name[sp]]][,c("object","new.sample.id")]
      }else{
        all.data <- rbind(all.data, ddf.dat.working[[species.name[sp]]][,c("object","new.sample.id")])
      }
    }
    #Update sample id's in obs.table
    obs.table <- merge(obs.table, all.data, by="object")[,c("object","Region.Label", "new.sample.id")]
    obs.table <- unique(obs.table)
    names(obs.table)[3] <- "Sample.Label"
    sample.table <- new.sample.table[,c("new.Sample.Label","Region.Label","Effort")]
    names(sample.table)[1] <- "Sample.Label"
  #RESAMPLE OBSERVATIONS WITHING SAMPLES   
  }else if(resample == "observations"){
    #for every ddf.dat.master element
    for(dat in seq(along = ddf.dat.master)){
      #get unique species codes
      species.code <- unique(ddf.dat.master[[dat]][[subset.variable]])
      for(sp in seq(along = species.code)){
        #select only for one species code
        temp.dat <- ddf.dat.master[[dat]][ddf.dat.master[[dat]][[subset.variable]] == species.code[sp],]
        #get sample ID's for data
        temp.dat <- merge(temp.dat, obs.table, by = 'object')
        sample.id <- unique(temp.dat$Sample.Label)
        #for every sample
        for(samp in seq(along = sample.id)){
          #get observation IDs
          object.id <- unique(temp.dat$object[which(temp.dat$Sample.Label == sample.id[samp])])
          #resample
          if(length(object.id) > 1){
            new.obs <- data.frame(object = sample(object.id, length(object.id), replace = TRUE))
          }else{  #sample() does something funny if it is only given one value or zero values 
            new.obs <- data.frame(object = object.id)
          }
          #find observations
          if(samp == 1){
            resample.result <- merge(new.obs, ddf.dat.master[[dat]], by = 'object')
          }else{    
            resample.result <- rbind(resample.result, merge(new.obs, ddf.dat.master[[dat]], by = 'object'))
          }         
        }#next sample
        if(sp == 1){
          ddf.dat.working[[species.name[dat]]] <- resample.result
        }else{
          ddf.dat.working[[species.name[dat]]] <- rbind(ddf.dat.working[[species.name[dat]]], resample.result)
        }       
        #ddf.dat.working[[dat]]$old.object.id <- ddf.dat.working[[species.name[dat]]]$object
        #ddf.dat.working[[species.name[dat]]]$object <- 1:nrow(ddf.dat.working[[species.name[dat]]])
      }#next species      
      renumber.duplicate.results <- renumber.duplicates(ddf.dat.working[[dat]], obs.table, double.observer)
      ddf.dat.working[[dat]] <- renumber.duplicate.results$ddf.dat
      obs.table <- renumber.duplicate.results$obs.table
    }#next dataset  
  }#end elseif 
  
  return(list(ddf.dat.working = ddf.dat.working, obs.table = obs.table, sample.table = sample.table))
}