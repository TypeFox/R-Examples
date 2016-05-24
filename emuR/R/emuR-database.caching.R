## Update cache of emuDB
## 
## Updates sqlite cache of loaded emuDB. This can be used
## to update changes to precached/loaded DBs as it only updates the deltas 
## in the cache which is considerably faster than reloading and therefore 
## recacheing the entire DB. This function is now called by load_emuDB if 
## load_emuDB finds a preexisting cache.
## @param emuDBhandle 
## @param verbose display infos
update_cache <- function(emuDBhandle, verbose = TRUE){

  DBconfig = load_DBconfig(emuDBhandle)
  
  # list sessions & bundles
  sessions = list_sessions(emuDBhandle)
  bundles = list_bundles(emuDBhandle)
  notUpdatedSessionDBI = list_sessionsDBI(emuDBhandle)
  notUpdatedBundlesDBI = list_bundlesDBI(emuDBhandle)
  
  # add column to sessions to track if already stored
  if(nrow(sessions) ==0){
    return()
  }
  sessions$stored = F
  
  progress = 0
  
  if(verbose){
    cat("INFO: Checking if cache needs update for ", nrow(bundles), " bundles...\n")
    pb = utils::txtProgressBar(min = 0, max = nrow(bundles), initial = progress, style=3)
    utils::setTxtProgressBar(pb, progress)
  }
  for(bndlIdx in 1:nrow(bundles)){
    sessionsDBI = list_sessionsDBI(emuDBhandle)
    
    bndl = bundles[bndlIdx,]
    
    # check if session needs to be added
    if(!bndl$session %in% sessionsDBI$name){
      add_sessionDBI(emuDBhandle, bndl$session)
    }
    
    # construct path to annotJSON
    annotFilePath = normalizePath(file.path(emuDBhandle$basePath, paste0(bndl$session, session.suffix), 
                                            paste0(bndl$name, bundle.dir.suffix), 
                                            paste0(bndl$name, bundle.annotation.suffix, '.json')))
    
    # calculate new MD5 sum of bundle annotJSON
    newMD5annotJSON = tools::md5sum(annotFilePath)
    names(newMD5annotJSON) = NULL
    
    # get old MD5 sum (NOTE: this returns an empty string if the bundle isn't present)
    oldMD5annotJSON = get_MD5annotJsonDBI(emuDBhandle, bndl$session, bndl$name)
    if(newMD5annotJSON != oldMD5annotJSON){
      # read annotJSON as charac 
      annotJSONchar = readChar(annotFilePath, file.info(annotFilePath)$size)
      
      # convert to bundleAnnotDFs
      bundleAnnotDFs = annotJSONcharToBundleAnnotDFs(annotJSONchar)
      
      if(length(oldMD5annotJSON) == 0){
        # add to bundle table
        add_bundleDBI(emuDBhandle, bndl$session, bndl$name, bundleAnnotDFs$annotates, bundleAnnotDFs$sampleRate, newMD5annotJSON)
      }else{
        # update bundle entry by
        # removing old bundle entry
        remove_bundleDBI(emuDBhandle, bndl$session, bndl$name)
        # and adding to bundle table
        add_bundleDBI(emuDBhandle, bndl$session, bndl$name, bundleAnnotDFs$annotates, bundleAnnotDFs$sampleRate, newMD5annotJSON)
        # and remove bundleAnnotDBI
        remove_bundleAnnotDBI(emuDBhandle, bndl$session, bndl$name)
      }
      # add to items, links, labels tables
      store_bundleAnnotDFsDBI(emuDBhandle, bundleAnnotDFs, bndl$session, bndl$name)
      
      # build redundat links and calc positions
      build_allRedundantLinks(emuDBhandle, bndl$session, bndl$name)
      calculate_postionsOfLinks(emuDBhandle)
    }
    
    # increase progress bar  
    progress=progress+1L
    if(verbose){
      utils::setTxtProgressBar(pb,progress)
    }
  }

  # remove superfluous sessions from session table
  superfluousSessions = dplyr::anti_join(notUpdatedSessionDBI, sessions, by = "name")
  if(nrow(superfluousSessions) > 0){
    for(sesIdx in 1:nrow(superfluousSessions)){
      remove_sessionDBI(emuDBhandle, superfluousSessions[sesIdx,])
    }
  }
  # remove superfluous bundles from bundle table and bundleAnnotDBI values from items, labels and links tables
  superfluousBundles = dplyr::anti_join(notUpdatedBundlesDBI, bundles, by = c("session", "name"))
  if(nrow(superfluousBundles) > 0){
    for(bndlIdx in 1:nrow(superfluousBundles)){
      remove_bundleDBI(emuDBhandle, superfluousBundles[bndlIdx,]$session, superfluousBundles[bndlIdx,]$name)
      remove_bundleAnnotDBI(emuDBhandle, sessionName = superfluousBundles[bndlIdx,]$session, bundleName = superfluousBundles[bndlIdx,]$name)
    }
  }
  
}

# FOR DEVELOPMENT 
# library('testthat') 
# test_file('tests/testthat/test_aaa_initData.R')
# test_file('tests/testthat/test_emuR-database.caching.R')
