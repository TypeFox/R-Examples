## Validates the DBI representation of bundle 
##
validate_bundleDBI <- function(emuDBhandle, session, bundle){
  
  DBconfig = load_DBconfig(emuDBhandle)
  
  # check that levels with same name are present
  levelNames <- DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT DISTINCT level FROM items WHERE db_uuid='", emuDBhandle$UUID, "'AND session = '", session,"' ",
                                                        "AND bundle ='", bundle, "'"))$level
  
  levelDefNames = sapply(DBconfig$levelDefinitions, function(l) l$name)
  delta1 = setdiff(levelNames, levelDefNames)
  delta2 = setdiff(levelDefNames, levelNames)
  
  if(length(delta1) != 0 || length(delta2) != 0){
    if(length(delta1) != 0){
      return(list(type = 'ERROR',
                  message = paste('Following levels where found that do not match any levelDefinition:', paste(delta1), ';',
                                  'in bundle:', bundle)))
    }else{
      warning("No items for levelDefinition where found for level:'", paste(delta2), "';", "in bundle:'", bundle , "'")
    }
  }
  
  # check that levels have same types
  levelTypes <- DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT DISTINCT level, type FROM items WHERE db_uuid='", emuDBhandle$UUID, "' AND session = '", session,"' ",
                                                        "AND bundle ='", bundle, "'"))$type
  
  levelDefTypes = sapply(DBconfig$levelDefinitions, function(l) l$type)

  delta1 = setdiff(levelTypes, levelDefTypes)
  delta2 = setdiff(levelDefTypes, levelTypes)
  
  if(length(delta1) != 0 || length(delta2) != 0){
    return(list(type = 'ERROR',
                message = paste('Following level types differ from those defined:', paste(levelNames[levelTypes != levelDefTypes], collapse = ', '), ';',
                                'in bundle:', bundle)))
  }  
  
  # validate sequence and overlaps in items of type SEGMENTS
  tmp <- DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT DISTINCT * FROM items WHERE session = '", session,"' ", 
                                                 "AND bundle ='", bundle, "' ",
                                                 "AND type = 'SEGMENT'"))
  
  #TODO: VALIDATE: SEQUENCE + OVERLAPS / LINKS'
  
  
  
  return(list(type = 'SUCCESS', 
              message = ''))
  
}




## FOR DEVELOPMENT
# library('testthat')
# test_file('tests/testthat/test_validate.R')
