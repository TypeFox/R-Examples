requireNamespace("RSQLite", quietly = T)

## Create emuDB DBconfig object from a TextGrid file
## 
## @param tgPath path to TextGrid file
## @param dbName name of the database
## @param basePath project base path
## @param tierNames character vector containing names of tiers to extract and convert. If NULL (the default) all
## tiers are converted.
## @return object of class emuDB.schema.db
## @import stringr uuid wrassp RSQLite
## @keywords internal
## 
create_DBconfigFromTextGrid = function(tgPath, dbName, basePath, tierNames = NULL){
  
  ####################
  # check parameters

  if(is.null(tgPath)) {
    stop("Argument tgPath (path to TextGrid file) must not be NULL\n")
  }
  
  if(is.null(dbName)) {
    stop("Argument dbName (name of new DB) must not be NULL\n")
  }
  
  #
  ####################
  
  # parse TextGrid  
  tgAnnotDFs = TextGridToBundleAnnotDFs(tgPath, sampleRate = 2000, name = "tmpBundleName", annotates = "tmpBundleName.wav") # sampleRate/name/annotates don't matter!! -> hardcoded
  
  # remove unwanted levels
  if(!is.null(tierNames)){
    # filter items
    tgAnnotDFs$items = dplyr::filter_(tgAnnotDFs$items, ~(level %in% tierNames))
    # filter labels
    tgAnnotDFs$labels = dplyr::filter_(tgAnnotDFs$labels, ~(name %in% tierNames))
  }
  
  levels = dplyr::distinct_(tgAnnotDFs$items, "level", .keep_all = TRUE)
  
  # create level definitions
  levelDefinitions = list()
  
  # generate defaultLvlOrder
  defaultLvlOrder=list()
  levIdx = 1  

  for(lineIdx in 1:nrow(levels)){
    lev = levels[lineIdx,]
    if(lev$type == 'SEGMENT' || lev$type == 'EVENT'){
      defaultLvlOrder[[length(defaultLvlOrder)+1L]]=lev$level
    }else{
      stop('Found levelDefinition that is not of type SEGMENT|EVENT while parsing TextGrid...this should not occur!')
    }
    # add new leveDef.
    levelDefinitions[[levIdx]] = list(name = lev$level, 
                                      type = lev$type, 
                                      attributeDefinitions = list(list(name = lev$level, type = "STRING")))
    levIdx = levIdx + 1
  }
  
  
  # create signalCanvas config
  sc = list(order = c("OSCI","SPEC"), 
            assign = list(), 
            contourLims = list())
  
  # create perspective
  defPersp = list(name = 'default', 
                  signalCanvases = sc, 
                  levelCanvases = list(order = defaultLvlOrder), 
                  twoDimCanvases = list(order = list()))
  # create EMUwebAppConfig 
  waCfg = list(perspectives = list(defPersp),
               activeButtons = list(saveBundle = TRUE,
                                    showHierarchy = TRUE))
  
  
  
  # generate full schema list
  dbSchema = list(name = dbName,
                  UUID = uuid::UUIDgenerate(),
                  mediafileExtension = 'wav',
                  ssffTrackDefinitions = list(),
                  levelDefinitions = levelDefinitions,
                  linkDefinitions = list(),
                  EMUwebAppConfig = waCfg)
  
  
  return(dbSchema)
}

# FOR DEVELOPMENT
# library('testthat')
# test_file('tests/testthat/test_aaa_initData.R')
# test_file('tests/testthat/test_emuR-create_DBconfigFromTextGrid.R')
