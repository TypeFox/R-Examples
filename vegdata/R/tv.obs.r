# db <- file.path(new_folder, dbs)
"tv.obs" <- function(db, tv_home, ...) {
    if(missing(tv_home)) tv_home <- tv.home() else if(tv_home != tv.home()) warning(paste("Given Turboveg root directory:", tv_home, "differ from the global root directory given by getOption('tv_home'):", getOption('tv_home')))
    # Observations
    obs <- read.dbf(file.path(tv_home, 'Data', db[1],'tvabund.dbf'))
    names(obs) <- TCS.replace(names(obs))
    ### Combine multiple databases
    # Would be more efficient with e.g. rbindlist but I am not sure if it worth the new package dependency
    if(length(db)>1) {
      cat('More than 1 database, trying to combine.\n')
    refl.1 <- tv.refl(db = db[1])
      for(i in 2:length(db)) {
      	refl.i <- tv.refl(db = db[i])
      	if(refl.1 != refl.i) {
          cat(db[1], 'vs.', db[i])
          stop('You are using different taxonomic reference lists in your databases!')
      	}
      	obs.tmp <- read.dbf(file.path(tv_home, 'Data', db[i],'tvabund.dbf'))
        names(obs.tmp) <- TCS.replace(names(obs.tmp))
        if(any(unique(obs$RELEVE_NR) %in% unique(obs.tmp$RELEVE_NR))) {
        		print(db[i])
            stop('Overlap of releve numbers between the databases!')
        	}	
    	  if(any(!names(obs) %in%  names(obs.tmp) ) | any(!names(obs.tmp) %in% names(obs))) {
  	      miss1 <- setdiff(colnames(obs), colnames(obs.tmp))
  	      miss2 <- setdiff(colnames(obs.tmp), colnames(obs))
  	      obs[, miss2] <- NA
  	      obs.tmp[, miss1] <- NA
  	      obs <- rbind(obs, obs.tmp) 
        } else obs <- rbind(obs, obs.tmp)
      }
    }
    class(obs) <- c('tv.obs','data.frame')
    obs
}
