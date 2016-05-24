
#####################################################
# functions used to build various path combinations
# plus helper functions

get_levelNameForAttributeName <- function(emuDBhandle, attributeName){
  DBconfig = load_DBconfig(emuDBhandle)
  for(lvlD in DBconfig$levelDefinitions){
    aNames = character(0)
    for(ad in lvlD$attributeDefinitions){
      aNames = c(aNames, ad[['name']])
      if(attributeName %in% aNames){
        return(lvlD[['name']])
      }
    }
  }
  return(NULL)
}


get_allAttributeNames<-function(emuDBhandle){
  DBconfig = load_DBconfig(emuDBhandle)
  aNames=character(0)
  for(lvlD in DBconfig$levelDefinitions){
    for(ad in lvlD$attributeDefinitions){
      aNames=c(aNames,ad$name)
    }
    
  }
  return(aNames)
}


get_linkLevelChildrenNames<-function(schema, superlevelName){
  chNames = character(0)
  for(ld in schema[['linkDefinitions']]){
    if(ld[['superlevelName']] == superlevelName){
      chNames=c(chNames, ld[['sublevelName']])
    }
  }
  return(chNames)
}

expand_linkPath <- function(p){
  expPath = list()
  pLen = length(p)
  if(pLen == 1){
    return(list())
  }
  expPath[[length(expPath)+1L]] = p
  expPath = c(expPath, expand_linkPath(p[1:(pLen-1)]))
  return(expPath)
}

## build all hierarchy paths including partial paths
## @return list containing paths and subpaths
build_allHierarchyPaths <- function(schema){
  extLds = list()
  for(ld in schema[['levelDefinitions']]){
    lName = ld[['name']]
    pathes = build_sublevelPathes(schema, lName)
    for(p in pathes){
      extLds = c(extLds, expand_linkPath(p))
    }
  }
  return(unique(extLds))
}


build_sublevelPathes <- function(DBconfig, levelName){
  pathes = list()
  chNames = get_linkLevelChildrenNames(DBconfig, levelName)
  if(length(chNames) == 0){
    pathes[[length(pathes) + 1L]] = c(levelName)
  }else{
    for(chName in chNames){
      chPathes = build_sublevelPathes(DBconfig, chName)
      for(chPath in chPathes){
        pathes[[length(pathes)+1L]] = c(levelName,chPath)
      }
    }
  }
  return(pathes)
}


build_levelPathes <- function(emuDBhandle){
  DBconfig = load_DBconfig(emuDBhandle)
  pathes = list()
  chNames = character(0)
  for(l in DBconfig$levelDefinitions){
    lPathes = build_sublevelPathes(DBconfig, l[['name']])
    pathes = c(pathes, lPathes)
  }
  return(pathes)
}

# get all paths through hierarchy connecting two levels
get_hierPathsConnectingLevels <- function(emuDBhandle, levelName1, levelName2){

  allHierPaths = build_allHierarchyPaths(load_DBconfig(emuDBhandle))
  
  conHierPaths = list()
  
  
  for(p in allHierPaths){
    # assume levelName1 is above levelName2
    if(p[1] == levelName1 & p[length(p)] == levelName2){
      conHierPaths[[length(conHierPaths) + 1]] = p
    }
    # assume levelName2 is above levelName1
    if(p[1] == levelName2 & p[length(p)] == levelName1){
      conHierPaths[[length(conHierPaths) + 1]] = p
    }
    
  }
  
  return(conHierPaths)
}



# builds "extended" link definitions
# lists link definitionsfor every possible directed connection between levels
# returns list of character vectors 
# the first element of each character vector contains the super level name of the levelDefinition,
# the follwing elements contain all exetnded linked sub level names  
build_extLinkDefinitions <- function(emuDBhandle){
  lds = list()
  pathes = build_levelPathes(emuDBhandle)
  for(p in pathes){
    pLen = length(p)
    for(i in 1:pLen){
      ld = character(0)
      for(j in i:pLen){
        ld = c(ld,p[j])
      }
      lds[[length(lds)+1L]] = ld
    }
  }
  return(lds)
}


find_segmentLevels<-function(emuDBhandle, attrName){
  lvlNm = get_levelNameForAttributeName(emuDBhandle, attrName)
  extLnkDefs = build_extLinkDefinitions(emuDBhandle)
  segLvlList=character(0)
  for(extLnkDef in extLnkDefs){
    if(extLnkDef[1]==lvlNm){
      for(trgLvlNm in extLnkDef[2:length(extLnkDef)]){
        
        trgLd=get_levelDefinition(emuDBhandle, trgLvlNm)
        if(trgLd['type']=='SEGMENT'){
          segLvlList=unique(c(segLvlList,trgLvlNm))
        }
      }
    }
  }
  return(segLvlList)
}

get_levelDefinition <- function(emuDBhandle, name){
  DBconfig = load_DBconfig(emuDBhandle)
  res = NULL
  for(ld in DBconfig$levelDefinitions){
    if(ld$name == name){
      res = ld
      break
    }
  }
  return(res)
}

###########################################
# DBconfig file handeling functions

## load function for _DBconfig.json file of emuDB
load_DBconfig <- function(emuDBhandle){
  dbCfgPath = file.path(emuDBhandle$basePath, paste0(emuDBhandle$dbName, database.schema.suffix))
  if(file.exists(dbCfgPath)){
    DBconfig = jsonlite::fromJSON(dbCfgPath, simplifyVector=FALSE)
  }else{
    stop(dbCfgPath, " does not seem to exist. This could be due to a bad 'name' entry in the DBconfig file. This field has to be the same as the name of the emuDB (directory & _DBconfig.json)")
  }
  return(DBconfig)
}

# store function for dbConfig
store_DBconfig <- function(emuDBhandle, dbConfig, basePath = NULL){
  if(is.null(basePath)){
    basePath = emuDBhandle$basePath
  }
  dbCfgPath = file.path(basePath, paste0(emuDBhandle$dbName, database.schema.suffix))
  json = jsonlite::toJSON(dbConfig, auto_unbox = TRUE, force = TRUE, pretty = TRUE)
  writeLines(json, dbCfgPath)
}


################################################################
################# CRUD DBconfig functions ######################
################################################################



###########################################
# CRUD operation for levelDefinitions

##' Add / List / Remove level definition to / of / from emuDB
##' 
##' Add / List / Remove database operation functions for level definitions. 
##' A level is a more general term for what is often referred to as a "tier". 
##' It is more general in the sense that people usually 
##' expect tiers to contain time information. Levels 
##' can either contain time information if they are of the 
##' type "EVENT" or of the type "SEGMENT" but are timeless 
##' if they are of the type "ITEM". For more information 
##' on the structural elements of an emuDB see \code{vignette(emuDB)}.
##' Note that a level cannot be removed, if it contains instances of annotation items
##' or if it is linked to another level.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param name name of level definition
##' @param type type of level definition ("SEGMENT","EVENT","ITEM")
##' @keywords emuDB database schema Emu
##' @name AddListRemoveLevelDefinitions
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' # add level called "Phonetic2" to the ae emuDB
##' # that could for example contain the transcriptions of a second annotator
##' add_levelDefinition(emuDBhandle = ae,
##'                     name = "Phonetic2",
##'                     type = "SEGMENT")
##'                     
##' # list level definition of ae emuDB
##' list_levelDefinitions(emuDBhandle = ae)
##' 
##' # remove newly added level definition
##' remove_levelDefinitions(emuDBhandle = ae,
##'                         name = "Phonetic2")
##' }
##' 
NULL

##' @rdname AddListRemoveLevelDefinitions
##' @export
add_levelDefinition<-function(emuDBhandle, name,
                              type){
  allowedTypes = c('ITEM', 'SEGMENT', 'EVENT')
  # precheck type 
  if(!(type %in% allowedTypes)){
    stop('Bad type given! Type has to be either ', paste(allowedTypes, collapse = ' | ') )
  }
  levelDefinition=list(name = name, type = type, 
                       attributeDefinitions = list(list(name = name, type = 'STRING')))
  dbConfig = load_DBconfig(emuDBhandle)
  # check if level definition (name) already exists 
  for(ld in dbConfig$levelDefinitions){
    if(ld$name == levelDefinition$name){
      stop("Level definition:", levelDefinition$name," already exists in database ", emuDBhandle$dbName)
    }
  }
  # add
  dbConfig$levelDefinitions[[length(dbConfig$levelDefinitions) + 1]] = levelDefinition
  
  store_DBconfig(emuDBhandle, dbConfig)
  
  invisible(NULL)
}


##' @rdname AddListRemoveLevelDefinitions
##' @export
list_levelDefinitions <- function(emuDBhandle){
  dbConfig = load_DBconfig(emuDBhandle)
  df <- data.frame(name = character(),
                   type = character(), 
                   nrOfAttrDefs = numeric(), 
                   stringsAsFactors = FALSE) 
  
  for(ld in dbConfig$levelDefinitions){
    df <- rbind(df, data.frame(name = ld$name, 
                               type = ld$type, 
                               nrOfAttrDefs = length(ld$attributeDefinitions),
                               stringsAsFactors = FALSE))
  }
  return(df)
}


##' @rdname AddListRemoveLevelDefinitions
##' @export
remove_levelDefinition<-function(emuDBhandle, name){
  
  dbConfig = load_DBconfig(emuDBhandle)
  # check if level definition (name)exists 
  if(!any(sapply(dbConfig$levelDefinitions, function(ld) ld[['name']] == name))){
    stop("Level definition:", name, " does not exist in database ", dbConfig$name)
  }
  # check if level is referenced by link defintion
  for(lkd in dbConfig$linkDefinitions){
    if(lkd[['superlevelName']] == name | lkd[['sublevelName']] == name){
      lkdStr = toString(lkd)
      stop("Cannot remove level definition ", name, ". It is referenced by link definition: ", lkdStr)
    }
  }
  
  # check if level is empty
  itemsDf = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT * FROM items i WHERE \
                        i.db_uuid='", emuDBhandle$UUID, "' AND i.level='", name, "'"))
  itemsCnt = nrow(itemsDf)
  if(itemsCnt > 0){
    stop("Level is not empty. Remove items first to delete level ", name)
  }
  
  # do removal
  newLvlDefs = list()
  for(lvlDef in dbConfig$levelDefinitions){
    if(lvlDef[['name']] != name){
      newLvlDefs[[length(newLvlDefs) + 1]] = lvlDef
    }
  }
  dbConfig$levelDefinitions = newLvlDefs
  
  store_DBconfig(emuDBhandle, dbConfig)
  
  return(invisible(NULL))
}

###################################################
# CRUD operations for attributeDefinitions

##' Add / List / Remove attribute definition to / of / from emuDB
##' 
##' Add / List / Remove database operation functions for attribute definition 
##' to / of / from an existing level definition
##' of a emuDB. Attribute
##' definitions can be viewed as definitions of
##' parallel labels for the annotational units (ITEMs) of the emuDB. 
##' Each level definition is required to have at least one 
##' default attribute definition that has the same name as the level definition
##' (automatically created by \code{\link{add_levelDefinition}}). For more 
##' information on the structural elements of an emuDB see \code{vignette(emuDB)}.
##' Note that as with level definitions, an attribute definition to a level cannot be removed,
##' if it contains labels in the emuDB.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param levelName name of level
##' @param name name of attributeDefinition
##' @param type type of attributeDefinition (currently only "STRING")
##' @keywords emuDB database DBconfig Emu 
##' @name AddListRemoveAttributeDefinitions
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' # add additional attribute definition to the "Phonetic" level
##' # of the ae emuDB that will contain the UTF8 IPA
##' # symbols of the phonetic transcriptions
##' add_attributeDefinition(emuDBhandle = ae,
##'                         levelName = "Phonetic",
##'                         name = "IPA-UTF8")
##'                         
##' # list attribute definitions for level "Word"
##' # of the ae emuDB
##' list_attributeDefinitions(emuDBhandle = ae, 
##'                           levelName = "Word")
##' 
##' # remove newly added attributeDefinition
##' remove_attributeDefinition(emuDBhandle = ae,
##'                            levelName = "Phonetic",
##'                            name = "IPA-UTF8")
##' }
##' 
NULL

##' @rdname AddListRemoveAttributeDefinitions
##' @export
add_attributeDefinition <- function(emuDBhandle, levelName, 
                                    name, type = "STRING"){
  if(type != "STRING"){
    stop("Currently only attributeDefinition of type 'STRING' allowed")
  }
  
  dbConfig = load_DBconfig(emuDBhandle)
  
  df = list_attributeDefinitions(emuDBhandle, levelName)
  
  
  if(!(name %in% df$name)){
    for(i in 1:length(dbConfig$levelDefinitions)){
      if(dbConfig$levelDefinitions[[i]]$name == levelName){
        dbConfig$levelDefinitions[[i]]$attributeDefinitions[[length(dbConfig$levelDefinitions[[i]]$attributeDefinitions) + 1]] = list(name = name, type = type)
        break
      }
    }
  }else{
    stop(paste0("attributeDefinition with name '", name, "' already present in level '", levelName, "'"))
  }
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
  
}


##' @rdname AddListRemoveAttributeDefinitions
##' @export
list_attributeDefinitions <- function(emuDBhandle, levelName){

  ld = get_levelDefinition(emuDBhandle, levelName)
  
  if(length(ld$attributeDefinitions) > 1){
    df = data.frame(name = character(), 
                    type = character(), 
                    hasLabelGroups = logical(), 
                    hasLegalLabels = logical(), 
                    stringsAsFactors = F)
    for(ad in ld$attributeDefinitions){
      df = rbind(df, df = data.frame(name = ad$name, 
                                     type = ad$type, 
                                     hasLabelGroups = !is.null(ad$labelGroups),
                                     hasLegalLabels = !is.null(ad$legalLabels),
                                     stringsAsFactors = F))
    }
  }else{
    df <- data.frame(name=ld$attributeDefinitions[[1]]$name, 
                     type=ld$attributeDefinitions[[1]]$type,
                     hasLabelGroups = !is.null(ld$attributeDefinitions[[1]]$labelGroups),
                     hasLegalLabels = !is.null(ld$attributeDefinitions[[1]]$legalLabels),
                     stringsAsFactors = F)
  }
  rownames(df) <- NULL
  return(df)
}


##' @rdname AddListRemoveAttributeDefinitions
##' @export
remove_attributeDefinition <- function(emuDBhandle, 
                                       levelName, 
                                       name){
  
  if(levelName == name){
    stop("Can not remove primary attributeDefinition (attributeDefinition with same name as level)")
  }
  
  dbConfig = load_DBconfig(emuDBhandle)
  
  ld = get_levelDefinition(emuDBhandle, levelName)
  
  # check if instances are present
  qRes = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT * FROM items AS it, labels AS lb WHERE ",
                                                   "it.db_uuid = lb.db_uuid AND ", 
                                                   "it.session = lb.session AND ", 
                                                   "it.bundle = lb.bundle AND ",
                                                   "it.item_id = lb.item_id AND ",
                                                   "it.level = '", levelName, "' AND ",
                                                   "lb.name = '", name, "'"))
  
  if(nrow(qRes) > 0){
    stop("Can not remove attributeDefinition if there are labels present")
  }else{
    levDefIdx = NULL
    for(i in 1:length(dbConfig$levelDefinitions)){
      if(dbConfig$levelDefinitions[[i]]$name == levelName){
        levDefIdx = i
        break
      }
    }
    
    for(i in 1:length(dbConfig$levelDefinitions[[levDefIdx]]$attributeDefinitions)){
      if(dbConfig$levelDefinitions[[levDefIdx]]$attributeDefinitions[[i]]$name == name){
        dbConfig$levelDefinitions[[levDefIdx]]$attributeDefinitions[[i]] = NULL
        break
      }
    }  
  }
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
  
}

###################################################
# CRUD operations for legalLabels

##' Set / Get / Remove legal labels of attributeDefinition of emuDB
##' 
##' Set / Get / Remove legal labels of a specific attributeDefinition of a emuDB. 
##' The legal labels are a character vector of strings
##' that specifies the labels that are legal (i.e. allowed / valid) for the given attribute. 
##' As the EMU-webApp won't allow the annotator to enter any labels that are not 
##' specified in this array, this is a simple way of assuring that a level 
##' has a consistent label set. For more information 
##' on the structural elements of an emuDB see \code{vignette(emuDB)}.
##' Note that defining legal labels for an attributeDefinition does not imply that the 
##' existing labels are checked for being 'legal' in the emuDB.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param levelName name of level
##' @param attributeDefinitionName name of attributeDefinition (can be and often is the level name)
##' @param legalLabels character vector of labels
##' @keywords emuDB database schema Emu
##' @name SetGetRemoveLegalLabels
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' legalPhoneticLabels = c("V", "m", "N", "s", "t", "H", "@:", "f", "r", 
##'                         "E", "n", "z", "S", "i:", "w", "@", "k", "I", "d", 
##'                         "db", "j", "u:", "dH", "l", "ai", "O", "D", "o:", "v")
##' 
##' # set legal labels of the 
##' # default "Phonetic" attributeDefinition of
##' # the "Phonetic" level of ae emuDB
##' set_legalLabels(emuDBhandle = ae, 
##'                 levelName = "Phonetic",
##'                 attributeDefinitionName = "Phonetic",
##'                 legalLabels = legalPhoneticLabels)
##' 
##' # get legal labels of the 
##' # default "Phonetic" attributeDefinition of
##' # the "Phonetic" level of ae emuDB
##' get_legalLabels(emuDBhandle = ae, 
##'                 levelName = "Phonetic", 
##'                 attributeDefinitionName = "Phonetic")
##'                 
##' 
##' # remove legal labels of the 
##' # default "Phonetic" attributeDefinition of
##' # the "Phonetic" level of ae emuDB
##' remove_legalLabels(emuDBhandle = ae, 
##'                    levelName = "Phonetic", 
##'                    attributeDefinitionName = "Phonetic")
##'                 
##' }
##' 
NULL

##' @rdname SetGetRemoveLegalLabels
##' @export
set_legalLabels <- function(emuDBhandle,
                            levelName,
                            attributeDefinitionName,
                            legalLabels){
  
  if(!is.null(legalLabels) & class(legalLabels) != "character"){
    stop("legalLables must be of class 'character'")
  }
  
  dbConfig = load_DBconfig(emuDBhandle)
  
  for(i in 1:length(dbConfig$levelDefinitions)){
    for(j in 1:length(dbConfig$levelDefinitions[[i]]$attributeDefinitions)){
      if(dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$name == attributeDefinitionName){
        dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$legalLabels = legalLabels
      }
    }
  }
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
  
}


##' @rdname SetGetRemoveLegalLabels
##' @export
get_legalLabels <- function(emuDBhandle,
                            levelName,
                            attributeDefinitionName){
  
  ld = get_levelDefinition(emuDBhandle, levelName)
  
  ll = NULL
  for(ad in ld$attributeDefinitions){
    if(ad$name == attributeDefinitionName){
      if(!is.null(ad$legalLabels)){
        ll = unlist(ad$legalLabels)
      }else{
        ll = NA
      }
    }
  }
  
  return(ll)
}


##' @rdname SetGetRemoveLegalLabels
##' @export
remove_legalLabels <- function(emuDBhandle,
                               levelName,
                               attributeDefinitionName){
  
  # remove by setting to NULL
  set_legalLabels(emuDBhandle,
                  levelName,
                  attributeDefinitionName,
                  legalLabels = NULL)
}

###################################################
# CRUD operations for attributeDefinition$labelGroups

##' Add / List / Remove labelGroup to / of / from attributeDefinition of emuDB
##' 
##' Add / List / Remove label group to / of / from a specific attribute definition. 
##' This label group can be used as a short hand  
##' to reference groups of labels specific
##' to an attribute definition (compared to global label groups that 
##' are added by \code{\link{add_labelGroup}}) in a 
##' \code{\link{query}}. A common example would be to
##' add a label group for something like the phonetic
##' category of nasals to be able reference them 
##' as "nasals" in a \code{\link{query}}. For more information 
##' on the structural elements of an emuDB see \code{vignette(emuDB)}.
##' 
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param levelName name of level
##' @param attributeDefinitionName name of attributeDefinition
##' @param labelGroupName name of label group
##' @param labelGroupValues character vector of labels
##' @keywords emuDB database schema Emu
##' @seealso add_labelGroup
##' @name AddListRemoveAttrDefLabelGroup
##' @examples
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' sampaNasals = c("m", "F", "n", "J", "N")
##' 
##' # add these values to the default Phonetic attribute
##' # definition of the Phonetic level of the ae emuDB
##' add_attrDefLabelGroup(emuDBhandle = ae,
##'                       levelName = "Phonetic",
##'                       attributeDefinitionName = "Phonetic",
##'                       labelGroupName = "sampaNasals",
##'                       labelGroupValues = sampaNasals)
##' 
##' # query the labelGroup
##' query(ae, "Phonetic=sampaNasals")
##' 
##' 
##' # list attribute definition label groups
##' # of attributeDefinition "Phonetic" of the level "Phonetic"
##' # of the ae emuDB
##' list_attrDefLabelGroups(emuDBhandle = ae, 
##'                         levelName = "Phonetic" , 
##'                         attributeDefinitionName = "Phonetic")
##' 
##' # remove the newly added attrDefLabelGroup
##' remove_attrDefLabelGroup(emuDBhandle = ae,
##'                          levelName = "Phonetic",
##'                          attributeDefinitionName = "Phonetic",
##'                          labelGroupName = "sampaNasals")
##' 
##' }
##' 
NULL

##' @rdname AddListRemoveAttrDefLabelGroup
##' @export
add_attrDefLabelGroup <- function(emuDBhandle,
                                  levelName,
                                  attributeDefinitionName, 
                                  labelGroupName,
                                  labelGroupValues){
  
  dbConfig = load_DBconfig(emuDBhandle)
  curLgs = list_attrDefLabelGroups(emuDBhandle, 
                                   levelName, 
                                   attributeDefinitionName)
  
  if(labelGroupName %in% curLgs$name){
    stop("labelGroupName '", labelGroupName ,"' already exists!")
  }
  for(i in 1:length(dbConfig$levelDefinitions)){
    for(j in 1:length(dbConfig$levelDefinitions[[i]]$attributeDefinitions)){
      if(dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$name == attributeDefinitionName){
        l = length(dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$labelGroups)
        dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$labelGroups[[l + 1]] = list(name = labelGroupName, 
                                                                                             values = labelGroupValues)
      }
    }
  }
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
}

##' @rdname AddListRemoveAttrDefLabelGroup
##' @export
list_attrDefLabelGroups <- function(emuDBhandle,
                                    levelName,
                                    attributeDefinitionName){

  ld = get_levelDefinition(emuDBhandle, levelName)
  
  df = data.frame(name = character(), 
                  values = character(),
                  stringsAsFactors = F)
  for(ad in ld$attributeDefinitions){
    if(ad$name == attributeDefinitionName){
      if(!is.null(ad$labelGroups)){
        for(lg in ad$labelGroups){
          df = rbind(df, data.frame(name = lg$name,
                                    values = paste0(lg$values, collapse = "; ") ))
        }
      }
    }
  }
  
  return(df)
}


##' @rdname AddListRemoveAttrDefLabelGroup
##' @export
remove_attrDefLabelGroup <- function(emuDBhandle,
                                     levelName,
                                     attributeDefinitionName, 
                                     labelGroupName){
  
  dbConfig = load_DBconfig(emuDBhandle)
  curLgs = list_attrDefLabelGroups(emuDBhandle, 
                                   levelName, 
                                   attributeDefinitionName)
  
  if(!labelGroupName %in% curLgs$name){
    stop("labelGroupName '", labelGroupName ,"' does not exists!")
  }
  
  for(i in 1:length(dbConfig$levelDefinitions)){
    for(j in 1:length(dbConfig$levelDefinitions[[i]]$attributeDefinitions)){
      if(dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$name == attributeDefinitionName){
        l = length(dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$labelGroups)
        dbConfig$levelDefinitions[[i]]$attributeDefinitions[[j]]$labelGroups[[l]] = NULL
      }
    }
  }
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
  
}

###################################################
# CRUD operations for linkDefinitions

##' Add / List / Remove linkDefinition to / of / from emuDB
##' 
##' Add / List / Remove new link definition to / of / from emuDB. A link definition
##' specifies the relationship between two levels, the
##' super-level and the sub-level. The entirety of all link 
##' definitions of a emuDB specifies the 
##' hierarchical structure of the database. For more information
##' on the structural elements of an emuDB see \code{vignette(emuDB)}.
##' 
##' Link type descriptions:
##' \itemize{
##' \item{\code{"ONE_TO_MANY"}}{A single ITEM of the super-level can be linked to multiple ITEMs of the sub-level}
##' \item{\code{"MANY_TO_MANY"}}{Multiple ITEMs of the super-level can be linked to multiple ITEMs of the sub-level}
##' \item{\code{"ONE_TO_ONE"}}{A single ITEM of the super-level can be linked to a single ITEM of the sub-level}
##' }
##' 
##' For all link types the rule applies that no links are allowed to cross any other links.
##' Further, a linkDefinition can not be removed, if there are links present in the emuDB.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param type type of linkDefinition (either \code{"ONE_TO_MANY"}, \code{"MANY_TO_MANY"} or \code{"ONE_TO_ONE"})
##' @param superlevelName name of super-level of linkDefinition
##' @param sublevelName name of sub-level of linkDefinition
##' @name AddListRemoveLinkDefinition
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded emuDB that was converted
##' # using the convert_TextGridCollection function called myTGcolDB
##' # (see ?load_emuDB and ?convert_TextGridCollection for more information)
##' 
##' # add link definition from super-level "Phoneme"
##' # to sub-level "Phonetic" of type "ONE_TO_MANY"
##' # for myTGcolDB emuDB
##' add_linkDefinition(emuDBhandle = myTGcolDB,
##'                    type = "ONE_TO_MANY",
##'                    superlevelName = "Phoneme",
##'                    sublevelName = "Phonetic")
##' 
##' # list link definitions for myTGcolDB emuDB
##' list_linkDefinitions(emuDBhandle = myTGcolDB)
##' 
##' # remove newly added link definition
##' remove_linkDefinition(emuDBhandle = myTGcolDB,
##'                       superlevelName = "Phoneme",
##'                       sublevelName = "Phonetic")
##' 
##' 
##' }
NULL

##' @rdname AddListRemoveLinkDefinition
##' @export
add_linkDefinition <- function(emuDBhandle, 
                               type,
                               superlevelName,
                               sublevelName){
  
  dbConfig = load_DBconfig(emuDBhandle)
  
  allowedTypes = c("ONE_TO_MANY", "MANY_TO_MANY", "ONE_TO_ONE")
  
  if(!type %in% allowedTypes){
    stop("Only the following types permitted: ", paste(allowedTypes, collapse = '; '))
  }
  
  curLds = list_linkDefinitions(emuDBhandle)
  
  # check if level is defined
  curLevs = list_levelDefinitions(emuDBhandle)
  if(!any(curLevs$name == superlevelName) | !any(curLevs$name == sublevelName)){
    stop("Either superlevelName or sublevelName are not defined")
  }
  
  
  # check if link between levels already exists
  if(any(curLds$superlevelName == superlevelName & curLds$sublevelName == sublevelName)){
    stop("linkDefinition already exists for superlevelName: '", 
         superlevelName, "' and sublevelName: '", sublevelName, "'")
  }
  
  l = length(dbConfig$linkDefinitions)
  dbConfig$linkDefinitions[[l + 1]] = list(type = type, 
                                           superlevelName = superlevelName,
                                           sublevelName = sublevelName)
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
  
}


##' @rdname AddListRemoveLinkDefinition
##' @export
list_linkDefinitions <- function(emuDBhandle){
  
  dbConfig = load_DBconfig(emuDBhandle)
  
  df = data.frame(type = character(),
                  superlevelName = character(),
                  sublevelName = character(),
                  stringsAsFactors = F)
  
  for(ld in dbConfig$linkDefinitions){
    df = rbind(df, data.frame(type = ld$type,
                              superlevelName = ld$superlevelName,
                              sublevelName = ld$sublevelName))
  }
  
  return(df)
  
}


##' @rdname AddListRemoveLinkDefinition
##' @export
remove_linkDefinition <- function(emuDBhandle, 
                                  superlevelName,
                                  sublevelName){
  
  dbConfig = load_DBconfig(emuDBhandle)
  curLds = list_linkDefinitions(emuDBhandle)
  
  # check if linkDef exists
  if(!any(curLds$superlevelName == superlevelName & curLds$sublevelName == sublevelName)){
    stop("No linkDefinition found for superlevelName '", superlevelName, 
         "' and sublevelName '", sublevelName, "'")
  }
  # check if links are present
  res = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT * FROM ",
                                                  "links ",
                                                  "INNER JOIN (SELECT * FROM items WHERE level = '", superlevelName, "' AND db_uuid = '", dbConfig$UUID, "') as superItems", 
                                                  "    ON links.from_id = superItems.item_id ",
                                                  "       AND links.db_uuid = superItems.db_uuid ",
                                                  "       AND links.session = superItems.session ",
                                                  "       AND links.bundle = superItems.bundle ",
                                                  "INNER JOIN (SELECT * FROM items WHERE level = '", sublevelName, "' AND db_uuid = '", dbConfig$UUID, "') as subItems", 
                                                  "    ON links.to_id = subItems.item_id ",
                                                  "       AND links.db_uuid = subItems.db_uuid ",
                                                  "       AND links.session = subItems.session ",
                                                  "       AND links.bundle = subItems.bundle ",
                                                  "WHERE links.db_uuid = '", emuDBhandle$UUID, "'"))
  
  if(nrow(res) != 0){
    stop("linkDefinition can not be remove as there are links present")
  }
  
  for(i in 1:length(dbConfig$linkDefinitions)){
    if(dbConfig$linkDefinitions[[i]]$superlevelName == superlevelName && dbConfig$linkDefinitions[[i]]$sublevelName == sublevelName){
      dbConfig$linkDefinitions[[i]] = NULL
      break
    }
  }
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
  
}

###################################################
# CRUD operations for ssffTrackDefinitions

##' Add / List / Remove ssffTrackDefinition to / from / of emuDB
##' 
##' Add / List / Remove ssffTrackDefinition to / from / of emuDB. 
##' An ssffTrack (often simply referred to as a track) references 
##' data that is stored in the Simple Signal File Format (SSFF) 
##' in the according bundle folders. The two most common types of data are:
##' \itemize{
##' \item{complementary data that was acquired during the recording 
##' such as data acquired during electromagnetic 
##' articulographic (EMA) or electropalatography (EPG) recordings;}
##' \item{derived data, i.e. data that was calculated from the original audio signal 
##' such as formant values and their bandwidths or the short-term Root Mean Square amplitude of the signal.}
##' }
##' For more information on the structural elements of an emuDB see \code{vignette(emuDB)}.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param name name of ssffTrackDefinition
##' @param columnName columnName of ssffTrackDefinition.
##' If the \code{onTheFlyFunctionName} parameter is set and columnName isn't, the
##' \code{columnName} will default to the first entry in \code{wrasspOutputInfos[[onTheFlyFunctionName]]$tracks}.
##' @param fileExtension fileExtension of ssffTrackDefinitions.
##' If the \code{onTheFlyFunctionName} parameter is set and fileExtension isn't, the
##' \code{fileExtension} will default to the first entry in \code{wrasspOutputInfos[[onTheFlyFunctionName]]$ext}.
##' @param onTheFlyFunctionName name of wrassp function to do on-the-fly calculation. If set to the name of a wrassp 
##' signal processing function, not only the emuDB schema is extended by the ssffTrackDefintion but also 
##' the track itself is calculated from the signal file and stored in the emuDB. See \code{names(wrasspOutputInfos)}
##' for a list of all the signal processing functions provided by the wrassp package.
##' @param onTheFlyParams a list of parameters that will be given to the function 
##' passed in by the onTheFlyFunctionName parameter. This list can easily be 
##' generated using the \code{\link{formals}} function on the according signal processing function 
##' provided by the wrassp package and then setting the
##' parameter one wishes to change.
##' @param onTheFlyOptLogFilePath path to optional log file for on-the-fly function
##' @param deleteFiles delete files that belong to ssffTrackDefinition on removal
##' @param verbose Show progress bars and further information
##' @param interactive ask user for confirmation
##' @seealso wrasspOutputInfos
##' @name AddListRemoveSsffTrackDefinition
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' # add ssff track definition to ae emuDB
##' # calculating the according SSFF files (.zcr) on-the-fly
##' # using the wrassp function "zcrana" (zero-crossing-rate analysis)
##' add_ssffTrackDefinition(emuDBhandle = ae,
##'                         name = "ZCRtrack",
##'                         onTheFlyFunctionName = "zcrana")
##'                         
##' # add ssff track definition to ae emuDB
##' # for SSFF files that will be added later (either
##' # by adding files to the emuDB using 
##' # the add_files() function or by calculating
##' # them using the according function provided 
##' # by the wrassp package)
##' add_ssffTrackDefinition(emuDBhandle = ae,
##'                         name = "formants",
##'                         columnName = "fm",
##'                         fileExtension = "fms")
##' 
##' # list ssff track definitions for ae emuDB
##' list_ssffTrackDefinitions(emuDBhandle = ae)
##' 
##' # remove newly added ssff track definition (does not delete 
##' # the actual .zrc files)
##' remove_ssffTrackDefinition <- function(emuDBhandle = ae, 
##'                                        name = "ZCRtrack")
##' 
##' }
##' 
NULL

##' @rdname AddListRemoveSsffTrackDefinition
##' @export
add_ssffTrackDefinition <- function(emuDBhandle, name,
                                    columnName = NULL, fileExtension = NULL, 
                                    onTheFlyFunctionName = NULL, onTheFlyParams = NULL, 
                                    onTheFlyOptLogFilePath = NULL,
                                    verbose = TRUE, interactive = TRUE){
  
  dbConfig = load_DBconfig(emuDBhandle)
  
  #########################
  # parameter checks
  
  # set columnName to fist tracks entry in wrasspOutputInfos if columnName is not set
  if(!is.null(onTheFlyFunctionName) && is.null(columnName)){
    columnName = wrasspOutputInfos[[onTheFlyFunctionName]]$tracks[1]
  }
  
  # set fileExtension to fist ext entry in wrasspOutputInfos if fileExtension is not set
  if(!is.null(onTheFlyFunctionName) && is.null(fileExtension)){
    fileExtension = wrasspOutputInfos[[onTheFlyFunctionName]]$ext[1]
  }
  
  
  # check if three main parameters are not null
  if(is.null(name) || is.null(columnName) || is.null(fileExtension)){
    stop('name, columnName, fileExtension have to be set!')
  }
  
  # check if onTheFlyFunctionName is set if onTheFlyParams is
  if(is.null(onTheFlyFunctionName) && !is.null(onTheFlyParams)){
    stop('onTheFlyFunctionName has to be set if onTheFlyParams is set!')
  }
  
  # check if both onTheFlyFunctionName and onTheFlyParams are set if onTheFlyOptLogFilePath is 
  if( !is.null(onTheFlyOptLogFilePath) && (is.null(onTheFlyFunctionName) || is.null(onTheFlyParams))){
    stop('Both onTheFlyFunctionName and onTheFlyParams have to be set for you to be able to use the onTheFlyOptLogFilePath parameter!')
  }
  
  curDefs = list_ssffTrackDefinitions(emuDBhandle)
  
  if(sum(curDefs$name == name) != 0){
    stop("ssffTrackDefinitions with name ", name ," already exists for emuDB: ", emuDBhandle$dbName, "!")
  }
  
  # calculate new files
  if(!is.null(onTheFlyFunctionName)){
    # check if files exist
    filesDf = list_files(emuDBhandle, fileExtension)
    ans = 'y'
    if(nrow(filesDf) != 0){
      fp = paste(emuDBhandle$basePath, paste0(fp$session, session.suffix), paste0(fp$bundle, bundle.dir.suffix), fp$file, sep = .Platform$file.sep)
      if(interactive){
        ans = readline(paste0("There are files present in '",emuDBhandle$dbName,"' that have the file extention '", 
                              fileExtension, "' Continuing will overwrite these files! Do you wish to proceed? (y/n) "))
      }
    }else{
      if(ans == 'y'){
        
        ###############################
        # set up function formals
        funcFormals = formals(onTheFlyFunctionName)
        funcFormals[names(onTheFlyParams)] = onTheFlyParams
        funcFormals$optLogFilePath = onTheFlyOptLogFilePath
        fp = list_files(emuDBhandle, dbConfig$mediafileExtension)
        funcFormals$listOfFiles = paste(emuDBhandle$basePath, paste0(fp$session, session.suffix), paste0(fp$bundle, bundle.dir.suffix), fp$file, sep = .Platform$file.sep)
        funcFormals$explicitExt = fileExtension
        
        # check if columnName is valid track
        if(!(columnName %in% wrasspOutputInfos[[onTheFlyFunctionName]]$tracks)){
          stop("'", columnName ,"' is not a column produced by '", onTheFlyFunctionName, "'! Please check wrasspOutputInfos for information on the tracks of each wrassp function.")
        }
        
        do.call(onTheFlyFunctionName, funcFormals)
      }else{
        stop('Aborted by user...')
      }
    }
  }
  
  # add new ssffTrackDefinition
  dbConfig$ssffTrackDefinitions[[length(dbConfig$ssffTrackDefinitions) + 1]] = list(name = name, 
                                                                                    columnName = columnName,
                                                                                    fileExtension = fileExtension)
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
}

##' @rdname AddListRemoveSsffTrackDefinition
##' @export
list_ssffTrackDefinitions <- function(emuDBhandle){
  dbConfig = load_DBconfig(emuDBhandle)
  df <- do.call(rbind, lapply(dbConfig$ssffTrackDefinitions, data.frame, stringsAsFactors=FALSE))
  return(df)
}


##' @rdname AddListRemoveSsffTrackDefinition
##' @export
remove_ssffTrackDefinition <- function(emuDBhandle, name, 
                                       deleteFiles = FALSE){
  
  dbConfig = load_DBconfig(emuDBhandle)
  
  # precheck if exists
  sDefs = list_ssffTrackDefinitions(emuDBhandle)  
  
  if(!(name %in% sDefs$name)){
    stop("No ssffTrackDefinitions found with name: '", name, "'")
  }
  # find end delete entry
  deletedDef = NULL
  for(i in 1:length(dbConfig$ssffTrackDefinitions)){
    if(dbConfig$ssffTrackDefinitions[[i]]$name == name){
      deletedDef = dbConfig$ssffTrackDefinitions[[i]]
      dbConfig$ssffTrackDefinitions[[i]] = NULL
      break
    }
  }
  
  # find and delete files
  if(deleteFiles){
    fp = list_files(emuDBhandle, deletedDef$fileExtension)
    fp = paste(emuDBhandle$basePath, paste0(fp$session, session.suffix), paste0(fp$bundle, bundle.dir.suffix), fp$file, sep = .Platform$file.sep)
    file.remove(fp)
  }
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
}

###################################################
# CRUD operations for global labelGroups

##' Add / List / Remove global labelGroup to / of / from emuDB
##' 
##' Add / List / Remove label group that can be used as a short hand  
##' to reference groups of labels that are globally defined
##' for the entire database (compared to attribute definition
##' specific label groups that 
##' are added by \code{\link{add_attrDefLabelGroup}}) in a 
##' \code{\link{query}}. A common example would be to
##' add a label group for something like the phonetic
##' category of nasals to be able to reference them 
##' as "nasals" in a \code{\link{query}}. 
##' In theory you could use a labelGroupName as a label instance within the
##' level, but since this could lead to serious confusion, it is better avoided.
##' For users transitioning from the legacy EMU system: Do not confuse a 
##' labelGroup with legal labels: a labelGroup 
##' had the unfortunate name 'legal labels' in the legacy EMU system.  
##' For more information on the structural elements of an emuDB 
##' see \code{vignette{emuDB}}.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param name name of label group
##' @param values character vector of labels
##' @keywords emuDB database schema Emu
##' @seealso add_attrDefLabelGroup
##' @name AddListRemoveLabelGroup
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' sampaNasals = c("m", "F", "n", "J", "N")
##' 
##' # add these values to the ae emuDB
##' # as a globally available labelGroup
##' add_labelGroup(emuDBhandle = ae,
##'                name = "sampaNasals",
##'                values = sampaNasals)
##' 
##' # query the labelGroup in the "Phonetic" level
##' query(emuDBhandle = ae, 
##'       query = "Phonetic == sampaNasals")
##' 
##' # query the labelGroup in the "Phoneme" level
##' query(emuDBhandle = ae, 
##'       query = "Phoneme == sampaNasals")
##' 
##' # list global label groups of ae emuDB
##' list_labelGroups(emuDBhandle = ae)
##' 
##' # remove the newly added labelGroup
##' remove_labelGroup(emuDBhandle = ae,
##'                   name = "sampaNasals")
##' }
##' 
NULL

##' @rdname AddListRemoveLabelGroup
##' @export
add_labelGroup <- function(emuDBhandle,
                           name,
                           values){
  
  dbConfig = load_DBconfig(emuDBhandle)
  curLgs = list_labelGroups(emuDBhandle)
  
  if(name %in% curLgs$name){
    stop("labelGroup with name '", name ,"' already exists!")
  }
  
  # add labelGroup
  dbConfig$labelGroups[[length(dbConfig$labelGroups) + 1]] = list(name = name, 
                                                                  values = values)
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
}


##' @rdname AddListRemoveLabelGroup
##' @export
list_labelGroups <- function(emuDBhandle){
  
  dbConfig = load_DBconfig(emuDBhandle)
  df = data.frame(name = character(),
                  values = character(),
                  stringsAsFactors = F)
  
  for(lg in dbConfig$labelGroups){
    df = rbind(df, data.frame(name = lg$name,
                              values = paste0(lg$values, collapse = "; ")))
  }
  
  return(df)
  
}


##' @rdname AddListRemoveLabelGroup
##' @export
remove_labelGroup <- function(emuDBhandle,
                              name){
  
  dbConfig = load_DBconfig(emuDBhandle)
  curLgs = list_labelGroups(emuDBhandle)
  
  if(!name %in% curLgs$name){
    stop("No labelGroup with name '", name ,"' found!")
  }
  
  for(i in 1:length(dbConfig$labelGroups)){
    if(dbConfig$labelGroups[[i]]$name == name){
      dbConfig$labelGroups[[i]] = NULL
    }
  }
  
  # store changes
  store_DBconfig(emuDBhandle, dbConfig)
}


# FOR DEVELOPMENT 
# library('testthat') 
# test_file('tests/testthat/test_aaa_initData.R')
# test_file('tests/testthat/test_emuR-database.DBconfig.R')