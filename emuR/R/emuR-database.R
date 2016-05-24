requireNamespace("stringr", quietly = T)
requireNamespace("uuid", quietly = T)
requireNamespace("wrassp", quietly = T)
requireNamespace("DBI", quietly = T)
requireNamespace("tidyjson", quietly = T)
requireNamespace("dplyr", quietly = T)

# constants

# API level of database object format
# increment this value if the internal database object format changes  
emuDB.apiLevel = 3L

# internalVars currently containing only server handle (should merge testingVars back into it as well)
.emuR_pkgEnv <- new.env()
assign("internalVars", list(testingVars = list(inMemoryCache = F)), envir = .emuR_pkgEnv)

#############################################
# file/folder suffixes of emuDB format

emuDB.suffix = '_emuDB'
session.suffix = '_ses'
bundle.dir.suffix = '_bndl'
bundle.annotation.suffix = '_annot'
database.schema.suffix = '_DBconfig.json'
database.cache.suffix = '_emuDBcache.sqlite'

#############################################
# create table / index definitions for DBI

database.DDL.emuDB = 'CREATE TABLE emu_db (
  uuid VARCHAR(36) NOT NULL,
  name TEXT,
  PRIMARY KEY (uuid)
);'

database.DDL.emuDB_session = 'CREATE TABLE session (
  db_uuid VARCHAR(36),
  name TEXT,
  PRIMARY KEY (db_uuid,name),
  FOREIGN KEY (db_uuid) REFERENCES emu_db(uuid) ON DELETE CASCADE
);'

database.DDL.emuDB_bundle = 'CREATE TABLE bundle (
  db_uuid VARCHAR(36),
  session TEXT,
  name TEXT,
  annotates TEXT,
  sample_rate FLOAT,
  md5_annot_json TEXT,
  PRIMARY KEY (db_uuid, session, name),
  FOREIGN KEY (db_uuid, session) REFERENCES session(db_uuid, name) ON DELETE CASCADE
);'

database.DDL.emuDB_items = 'CREATE TABLE items (
  db_uuid VARCHAR(36),
  session TEXT,
  bundle TEXT,
  item_id INTEGER,
  level TEXT,
  type TEXT,
  seq_idx INTEGER,
  sample_rate FLOAT,
  sample_point INTEGER,
  sample_start INTEGER,
  sample_dur INTEGER,
  PRIMARY KEY (db_uuid, session, bundle, item_id),
  FOREIGN KEY (db_uuid, session, bundle) REFERENCES bundle(db_uuid, session, name) ON DELETE CASCADE
);'

# database.DDL.emuDB_items_idx1 = 'CREATE INDEX IF NOT EXISTS items_idx1 ON items(item_id, bundle, session, db_uuid)'

# Important note:
# The primary key of items contains more columns then needed to identify a particular item.
# PRIMARY KEY (db_uuid,session,bundle,item_id) would be sufficient but the extended primary key 
# is necessary to speed up the build_redundnatLinksForPathes SQL query.
# It did not work to create an index like the one in the comment line below.
# It seems the query always uses the index of the primary key.
#database.DDL.emuDB_itemsIdx='CREATE UNIQUE INDEX items_level_idx ON items(db_uuid,session,bundle,level,item_id,type)'

database.DDL.emuDB_labels = 'CREATE TABLE labels (
  db_uuid VARCHAR(36),
  session TEXT,
  bundle TEXT,
  item_id INTEGER,
  label_idx INTEGER,
  name TEXT,
  label TEXT,
  PRIMARY KEY (db_uuid, session, bundle, item_id, label_idx),
  FOREIGN KEY (db_uuid, session, bundle) REFERENCES bundle(db_uuid, session, name) ON DELETE CASCADE
);'

# database.DDL.emuDB_labels_idx1 = 'CREATE INDEX IF NOT EXISTS labels_idx1 ON labels(item_id, bundle, session, db_uuid)'

database.DDL.emuDB_links = 'CREATE TABLE links (
  db_uuid VARCHAR(36) NOT NULL,
  session TEXT,
  bundle TEXT,
  from_id INTEGER,
  to_id INTEGER,
  label TEXT,
  FOREIGN KEY (db_uuid, session, bundle) REFERENCES bundle(db_uuid, session, name) ON DELETE CASCADE
);'

database.DDL.emuDB_links_both_ids_idx = 'CREATE INDEX IF NOT EXISTS links_both_ids_idx ON links(db_uuid, session, bundle, from_id, to_id)'
database.DDL.emuDB_links_to_id_idx = 'CREATE INDEX IF NOT EXISTS links_to_id_idx ON links(db_uuid, session, bundle, to_id)'

database.DDL.emuDB_linksTmp = 'CREATE TEMP TABLE links_tmp (
  db_uuid VARCHAR(36) NOT NULL,
  session TEXT,
  bundle TEXT,
  from_id INTEGER,
  to_id INTEGER,
  label TEXT
);'

database.DDL.emuDB_linksTmpIdx = 'CREATE INDEX IF NOT EXISTS links_tmp_idx ON links_tmp(db_uuid, session, bundle, from_id, to_id)'

database.DDL.emuDB_linksExt = 'CREATE TABLE links_ext (
  db_uuid VARCHAR(36) NOT NULL,
  session TEXT,
  bundle TEXT,
  from_id INTEGER,
  to_id INTEGER,
  seq_idx INTEGER,
  to_level TEXT,
  type TEXT,
  to_seq_idx INTEGER,
  to_seq_len INTEGER,
  label TEXT,
  FOREIGN KEY (db_uuid, session, bundle) REFERENCES bundle(db_uuid, session, name) ON DELETE CASCADE
);'


database.DDL.emuDB_linksExtIdx = 'CREATE INDEX IF NOT EXISTS links_ext_idx ON links_ext(db_uuid,session,bundle,from_id,to_id,to_level,type)'

# this should be a temp table
database.DDL.emuDB_linksExtTmp = 'CREATE TEMP TABLE links_ext_tmp (
  db_uuid VARCHAR(36) NOT NULL,
  session TEXT,
  bundle TEXT,
  from_id INTEGER,
  to_id INTEGER,
  seq_idx INTEGER,
  to_level TEXT,
  type TEXT,
  to_seq_idx INTEGER,
  to_seq_len INTEGER,
  label TEXT
);'

database.DDL.emuDB_linksExtTmpIdx = 'CREATE INDEX IF NOT EXISTS links_ext_tmp_idx ON links_ext_tmp(db_uuid,session,bundle,from_id,to_id,to_level,type)'

# this should be a temp table
database.DDL.emuDB_linksExtTmp2 = 'CREATE TEMP TABLE links_ext_tmp2 (
  db_uuid VARCHAR(36) NOT NULL,
  session TEXT,
  bundle TEXT,
  from_id INTEGER,
  to_id INTEGER,
  seq_idx INTEGER,
  to_level TEXT,
  type TEXT,
  to_seq_idx INTEGER,
  to_seq_len INTEGER,
  label TEXT
);'

database.DDL.emuDB_linksExtTmpIdx2 = 'CREATE INDEX IF NOT EXISTS links_ext_tmp2_idx ON links_ext_tmp2(db_uuid,session,bundle,from_id,to_id,to_level,type)'

####################################
######### DBI functions ############
####################################

####################################
# init functions (create tables and indices)

initialize_emuDbDBI <- function(emuDBhandle, createTables=TRUE, createIndices=TRUE){
  # check of old tables are present and rename them
  if(DBI::dbExistsTable(emuDBhandle$connection, "emuDB")){
    cat("INFO: Depricated cache tables found. Deleting these and recreating SQL cache that adheres to new DB schema definition.\n")
    allTableNames = DBI::dbListTables(emuDBhandle$connection)
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "items"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "labels"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "links"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "linksTmp"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "linksExt"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "linksExtTmp"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "linksExtTmp2"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "bundle"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "session"))
    DBI::dbGetQuery(emuDBhandle$connection, paste0("DROP TABLE IF EXISTS ", "emuDB"))
  }
    
  if(createTables & !DBI::dbExistsTable(emuDBhandle$connection, "emu_db")){
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_session)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_bundle)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_items)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_labels)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_links)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksExt)
  }
  if(createIndices){
    create_emuDBindicesDBI(emuDBhandle)
  }
}

create_emuDBindicesDBI<-function(emuDBhandle){
  DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_links_both_ids_idx)
  DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_links_to_id_idx)
  # DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_labels_idx1)
  # DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksIdx)
  DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksExtIdx)
}


####################################
# emuDB table DBI functions

add_emuDbDBI <- function(emuDBhandle){
  dbSqlInsert = paste0("INSERT INTO emu_db(uuid,name) VALUES('", emuDBhandle$UUID, "','", emuDBhandle$dbName, "')")
  DBI::dbGetQuery(emuDBhandle$connection, dbSqlInsert)
}

get_emuDbDBI <- function(emuDBhandle){
  query = paste0("SELECT * FROM emu_db WHERE uuid='", emuDBhandle$UUID, "'")
  res <- DBI::dbGetQuery(emuDBhandle$connection, query)
  return(res)
}


####################################
# session table DBI functions

add_sessionDBI <- function(emuDBhandle, sessionName){
  insertSessionSql = paste0("INSERT INTO session(db_uuid, name) VALUES('", emuDBhandle$UUID,"','", sessionName, "')")
  DBI::dbGetQuery(emuDBhandle$connection, insertSessionSql)
}

list_sessionsDBI <- function(emuDBhandle){
  dbs=DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT name FROM session WHERE db_uuid='", emuDBhandle$UUID, "'"))
  return(dbs)
}


remove_sessionDBI <- function(emuDBhandle, sessionName){
  DBI::dbGetQuery(emuDBhandle$connection, paste0("DELETE FROM session WHERE ", "db_uuid='", emuDBhandle$UUID, "' AND name='", sessionName, "'"))
}

####################################
# bundle table DBI functions

add_bundleDBI <- function(emuDBhandle, sessionName, name, annotates, sampleRate, MD5annotJSON){
  insertBundleSql = paste0("INSERT INTO bundle(db_uuid, session, name, annotates, sample_rate, md5_annot_json) VALUES('", 
                           emuDBhandle$UUID, "', '", sessionName, "', '", name, "', '", annotates, "', ", sampleRate, ", '", MD5annotJSON, "')")
  DBI::dbGetQuery(emuDBhandle$connection, insertBundleSql)
}

list_bundlesDBI <- function(emuDBhandle, sessionName = NULL){
  if(is.null(sessionName)){
    bundle = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT session, name FROM bundle WHERE db_uuid='", emuDBhandle$UUID, "'"))
  }else{
    bundle = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT session, name FROM bundle WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "'"))
  }
  return(bundle)
}

remove_bundleDBI <- function(emuDBhandle, sessionName, name){
  DBI::dbGetQuery(emuDBhandle$connection, paste0("DELETE FROM bundle WHERE ", "db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND name='", name, "'"))
}

# MD5annotJSON
get_MD5annotJsonDBI <- function(emuDBhandle, sessionName, name){
  MD5annotJSON = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT md5_annot_json as md5 FROM bundle WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND name='", name, "'"))$md5
  if(length(MD5annotJSON) == 0){
    MD5annotJSON = ""
  }
  return(MD5annotJSON)
}

####################################
# items, links, labels DBI functions

store_bundleAnnotDFsDBI <- function(emuDBhandle, bundleAnnotDFs, sessionName, 
                                    bundleName) {
  
  # insert items table entries (fist exanding it with db_uuid, session and bundle columns)
  if(nrow(bundleAnnotDFs$items) > 0){
    bundleAnnotDFs$items = data.frame(db_uuid = emuDBhandle$UUID, 
                                      session = sessionName,
                                      bundle = bundleName,
                                      bundleAnnotDFs$items)

    DBI::dbWriteTable(emuDBhandle$connection, "items", bundleAnnotDFs$items, append = T, row.names = F)
  }
  
  # insert labels table entries (fist exanding it with db_uuid, session and bundle columns)
  if(nrow(bundleAnnotDFs$labels) > 0){
    bundleAnnotDFs$labels =  data.frame(db_uuid = emuDBhandle$UUID, 
                                        session = sessionName,
                                        bundle = bundleName,
                                        bundleAnnotDFs$labels)
    
    DBI::dbWriteTable(emuDBhandle$connection, "labels", bundleAnnotDFs$labels, append = T, row.names = F)
  }
  
  # insert links table entries (fist exanding it with db_uuid, session and bundle columns)
  if(nrow(bundleAnnotDFs$links) > 0){
    bundleAnnotDFs$links =  data.frame(db_uuid = emuDBhandle$UUID,
                                       session = sessionName,
                                       bundle = bundleName,
                                       bundleAnnotDFs$links,
                                       label = NA)
    
    DBI::dbWriteTable(emuDBhandle$connection, "links", bundleAnnotDFs$links, append = T, row.names = F)
  }
}

load_bundleAnnotDFsDBI <- function(emuDBhandle, sessionName, bundleName){

  DBconfig = load_DBconfig(emuDBhandle)
  levelDefs = list_levelDefinitions(emuDBhandle)
  # meta infos
  annotates = paste0(bundleName, ".", DBconfig$mediafileExtension)
  sampleRateQuery = paste0("SELECT sample_rate FROM bundle WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND name='", bundleName,"'")
  sampleRate = DBI::dbGetQuery(emuDBhandle$connection, sampleRateQuery)$sample_rate
  
  # items
  itemsQuery = paste0("SELECT item_id, level, type, seq_idx, sample_rate, sample_point, sample_start, sample_dur FROM items WHERE db_uuid='", 
                      emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName,"' ORDER BY level, seq_idx")
  items = DBI::dbGetQuery(emuDBhandle$connection, itemsQuery)
  # reorder items to match DBconfig
  items = items[order(match(items$level,levelDefs$name)),]
  
  # labels 
  labelsQuery = paste0("SELECT item_id, label_idx, name, label FROM labels WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName,"'")
  labels = DBI::dbGetQuery(emuDBhandle$connection, labelsQuery)
  
  # links 
  
  linksQuery = paste0("SELECT from_id, to_id, label FROM links WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName,"'")
  links = DBI::dbGetQuery(emuDBhandle$connection, linksQuery)
  
  
  return(list(name = bundleName, annotates = annotates, sampleRate = sampleRate, items = items, links = links, labels = labels))
}


remove_bundleAnnotDBI<-function(emuDBhandle, sessionName, bundleName){
  cntSqlQuery=paste0("SELECT * FROM items WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName,"'")
  res<-DBI::dbGetQuery(emuDBhandle$connection, cntSqlQuery)
  
  delSqlQuery=paste0("DELETE FROM items WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName, "'")
  DBI::dbGetQuery(emuDBhandle$connection, delSqlQuery)
  
  delSqlQuery=paste0("DELETE FROM labels WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName,"'")
  DBI::dbGetQuery(emuDBhandle$connection, delSqlQuery)
  
  delSqlQuery=paste0("DELETE FROM links WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName, "'")
  DBI::dbGetQuery(emuDBhandle$connection, delSqlQuery)
  
  cntSqlQuery=paste0("SELECT * FROM links_ext WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName, "'")
  res<-DBI::dbGetQuery(emuDBhandle$connection, cntSqlQuery)
  
  delSqlQuery=paste0("DELETE FROM links_ext WHERE db_uuid='", emuDBhandle$UUID, "' AND session='", sessionName, "' AND bundle='", bundleName,"'")
  DBI::dbGetQuery(emuDBhandle$connection,delSqlQuery)
}



###################################################
# create redundant links functions

create_tmpTablesForBuildingRedLinks <- function(emuDBhandle){
  if(!"links_tmp" %in% DBI::dbListTables(emuDBhandle$connection)){
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksTmp)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksTmpIdx)
  }
  if(!"links_ext_tmp" %in% DBI::dbListTables(emuDBhandle$connection)){
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksExtTmp)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksExtTmpIdx)
    }
  if(!"links_ext_tmp2" %in% DBI::dbListTables(emuDBhandle$connection)){ 
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksExtTmp2)
    DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuDB_linksExtTmpIdx2)
    }
}

drop_tmpTablesForBuildingRedLinks <- function(emuDBhandle){
  if("links_tmp" %in% DBI::dbListTables(emuDBhandle$connection)) DBI::dbGetQuery(emuDBhandle$connection, "DROP TABLE links_tmp")
  if("links_ext_tmp" %in% DBI::dbListTables(emuDBhandle$connection)) DBI::dbGetQuery(emuDBhandle$connection, "DROP TABLE links_ext_tmp")
  if("links_ext_tmp2" %in% DBI::dbListTables(emuDBhandle$connection)) DBI::dbGetQuery(emuDBhandle$connection, "DROP TABLE links_ext_tmp2")
}

## Legacy EMU and query functions link collections contain links for each possible connection between levels
## We consider links that do not follow link definition constraints as redundant and therefore we remove them from the
## link data model
build_allRedundantLinks <- function(emuDBhandle, sessionName=NULL, bundleName=NULL){
  
  hierarchyPaths = build_allHierarchyPaths(load_DBconfig(emuDBhandle))
  
  return(build_redundantLinksForPaths(emuDBhandle, hierarchyPaths, sessionName, bundleName) )
}


build_redundantLinksForPaths <- function(emuDBhandle, hierarchyPaths, sessionName='0000', bundleName=NULL){
  
  # create tmp tables if not available
  create_tmpTablesForBuildingRedLinks(emuDBhandle)
  # delete any previous redundant links (just to be safe)
  DBI::dbGetQuery(emuDBhandle$connection, 'DELETE FROM links_tmp')
  
  hierarchyPathsLen = length(hierarchyPaths)
  if(hierarchyPathsLen > 0){
    
    sqlQuery = "INSERT INTO links_tmp(db_uuid,session,bundle,from_id,to_id,label) SELECT DISTINCT f.db_uuid,f.session,f.bundle,f.item_id AS from_id,t.item_id AS to_id, NULL AS label FROM items f,items t"
    sqlQuery = paste0(sqlQuery," WHERE f.db_uuid='", emuDBhandle$UUID, "' AND f.db_uuid=t.db_uuid AND f.session=t.session AND f.bundle=t.bundle AND ")
    #sqlQuery=paste0(sqlQuery," WHERE f.db_uuid=t.db_uuid AND f.bundle=t.bundle AND f.session=t.session AND ")
    
    if(!is.null(sessionName) & !is.null(bundleName)){
      # only for one bundle
      sqlQuery = paste0(sqlQuery,"f.session='",sessionName,"' AND f.bundle='",bundleName,"' AND ")
    }
    
    sqlQuery=paste0(sqlQuery,' (')
    
    # build query for each partial path
    for(i in 1:hierarchyPathsLen){
      hp = hierarchyPaths[[i]]
      #cat("Path: ",hp,"\n")
      hpLen = length(hp)
      sHp = hp[1]
      eHp = hp[hpLen]
      sqlQuery = paste0(sqlQuery, "(f.level='", sHp, "' AND t.level='", eHp, "'" )
      sqlQuery = paste0(sqlQuery, " AND EXISTS (SELECT l1.* FROM ")
      for(li in 1:(hpLen - 1)){
        sqlQuery = paste0(sqlQuery, 'links l', li)
        if(li < (hpLen - 1)){
          sqlQuery = paste0(sqlQuery, ',')
        }
      }
      if(hpLen > 2){
        for(ii in 2:(hpLen-1)){
          sqlQuery=paste0(sqlQuery,',items i',ii)
        }
      }
      sqlQuery=paste0(sqlQuery," WHERE ")
      if(hpLen==2){
        sqlQuery=paste0(sqlQuery,"l1.db_uuid=f.db_uuid AND l1.db_uuid=t.db_uuid AND l1.session=f.session AND l1.session=t.session AND l1.bundle=f.bundle AND l1.bundle=t.bundle AND f.item_id=l1.from_id AND t.item_id=l1.to_id")
        
      }else{
        # TODO start and end connection
        # from start to first in-between item 
        eHp=hp[2]
        sqlQuery=paste0(sqlQuery,"l1.db_uuid=f.db_uuid AND l1.db_uuid=i2.db_uuid AND l1.session=f.session AND l1.session=i2.session AND l1.bundle=f.bundle AND l1.bundle=i2.bundle AND f.item_id=l1.from_id AND i2.item_id=l1.to_id AND f.level='",sHp,"' AND i2.level='",eHp,"' AND ")
        if(hpLen>3){
          for(j in 2:(hpLen-2)){
            sHp=hp[j]
            eHp=hp[j+1L] 
            sqlQuery=paste0(sqlQuery,"l",j,".db_uuid=i",j,".db_uuid AND l",j,".db_uuid=i",(j+1),".db_uuid AND l",j,".session=i",j,".session AND l",j,".session=i",(j+1),".session AND l",j,".bundle=i",j,".bundle AND l",j,".bundle=i",(j+1),".bundle AND i",j,".item_id=l",j,".from_id AND i",(j+1L),".item_id=l",j,".to_id AND i",j,".level='",sHp,"' AND i",(j+1L),".level='",eHp,"' AND ")
          }
        }
        # from last in-between item to end item
        sHp=hp[(hpLen-1)]
        eHp=hp[hpLen]
        
        j=hpLen-1
        sqlQuery=paste0(sqlQuery,"l",j,".db_uuid=i",j,".db_uuid AND l",j,".db_uuid=t.db_uuid AND l",j,".session=i",j,".session AND l",j,".session=t.session AND l",j,".bundle=i",j,".bundle AND l",j,".bundle=t.bundle AND i",j,".item_id=l",j,".from_id AND t.item_id=l",j,".to_id AND i",j,".level='",sHp,"' AND t.level='",eHp,"'")
      }
      sqlQuery=paste0(sqlQuery,"))")
      if(i<hierarchyPathsLen){
        sqlQuery=paste0(sqlQuery," OR ")
      }
    }
    sqlQuery=paste0(sqlQuery,")")
    # since version 2.8.x of sqlite the query is very slow without indices
    DBI::dbGetQuery(emuDBhandle$connection, sqlQuery)
  }
}

calculate_postionsOfLinks<-function(emuDBhandle){
  
  # for all position related functions we need to calculate the sequence indices of dominated items grouped to one dominance item 
  # Extend links table with sequence index of the targeted (dominated) item
  DBI::dbGetQuery(emuDBhandle$connection,"DELETE FROM links_ext_tmp")
  
  DBI::dbGetQuery(emuDBhandle$connection,"INSERT INTO links_ext_tmp(db_uuid,session,bundle,from_id,to_id,seq_idx,to_level,type,label) SELECT k.db_uuid,k.session,k.bundle,k.from_id,k.to_id,i.seq_idx,i.level AS to_level,i.type,NULL AS label FROM links_tmp k,items i WHERE i.db_uuid=k.db_uuid AND i.session=k.session AND i.bundle=k.bundle AND k.to_id=i.item_id")
  
  # extend links table with relative sequence index
  DBI::dbGetQuery(emuDBhandle$connection,"INSERT INTO links_ext_tmp2(db_uuid,session,bundle,seq_idx,from_id,to_id,to_level,type,label,to_seq_idx) SELECT k.db_uuid,k.session,k.bundle,k.seq_idx,k.from_id,k.to_id,k.to_level,k.type,k.label,k.seq_idx-(SELECT MIN(m.seq_idx) FROM links_ext_tmp m WHERE m.from_id=k.from_id AND m.db_uuid=k.db_uuid AND m.session=k.session AND m.bundle=k.bundle AND k.to_level=m.to_level GROUP BY m.db_uuid,m.session,m.bundle,m.from_id,m.to_level) AS to_seq_idx FROM links_ext_tmp k")
  
  DBI::dbGetQuery(emuDBhandle$connection,"DELETE FROM links_ext_tmp")
  
  # Add length of dominance group sequence
  DBI::dbGetQuery(emuDBhandle$connection,"INSERT INTO links_ext(db_uuid,session,bundle,seq_idx,from_id,to_id,to_seq_idx,to_level,type,label,to_seq_len) SELECT k.db_uuid,k.session,k.bundle,k.seq_idx,k.from_id,k.to_id,k.to_seq_idx,k.to_level,k.type,k.label,(SELECT MAX(m.seq_idx)-MIN(m.seq_idx)+1 FROM links_ext_tmp2 m WHERE m.from_id=k.from_id AND m.db_uuid=k.db_uuid AND m.session=k.session AND m.bundle=k.bundle AND k.to_level=m.to_level GROUP BY m.db_uuid,m.session,m.bundle,m.from_id,m.to_level) AS to_seq_len FROM links_ext_tmp2 k")
  DBI::dbGetQuery(emuDBhandle$connection,"DELETE FROM links_ext_tmp2")
  
  # remove temporary tables
  drop_tmpTablesForBuildingRedLinks(emuDBhandle)
}

##########################################
################# emuDB ##################
##########################################

#############################################
# function that use emuDB files (vs. DBI)

##' List sessions of emuDB
##' @description List session names of emuDB
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @return data.frame object with session names
##' @export
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB
##' # (see ?load_emuDB for more information)
##' 
##' # list all sessions of ae emuDB
##' list_sessions(emuDBhandle = ae)
##' 
##' }
##' 
list_sessions <- function(emuDBhandle){
  sesPattern = paste0("^.*", session.suffix ,"$")
  sesDirs = dir(emuDBhandle$basePath, pattern = sesPattern)
  sesDirs = gsub(paste0(session.suffix, "$"), "", sesDirs)
  return(data.frame(name = sesDirs, stringsAsFactors = F))
}

##' List bundles of emuDB
##' 
##' List all bundles of emuDB or of particular session.
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param session optional session
##' @return data.frame object with columns session and name of bundles
##' @export
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB
##' # (see ?load_emuDB for more information)
##' 
##' # list bundles of session "0000" of ae emuDB
##' list_bundles(emuDBhandle = ae,
##'              session = "0000")
##' 
##' }
##' 
list_bundles <- function(emuDBhandle, session=NULL){
  sesDf = list_sessions(emuDBhandle)
  if(!is.null(session)){
    sesDf = dplyr::filter_(sesDf, ~(name == session))
  }
  bndlPattern = paste0("^.*", bundle.dir.suffix ,"$")
  res = data.frame(session = character(), name = character(), stringsAsFactors = F)
  
  
  for(ses in sesDf$name){
    bndlDirs = dir(file.path(emuDBhandle$basePath, paste0(ses, session.suffix)), pattern = bndlPattern)
    bndlNames = gsub(paste0(bundle.dir.suffix, "$"), "", bndlDirs)
    if(length(bndlNames) > 0){
      res = rbind(res, data.frame(session = ses, name = bndlNames, stringsAsFactors = F))
    }
  }
  return(res)
}



rewrite_allAnnots <- function(emuDBhandle, verbose=TRUE){
  
  bndls = list_bundles(emuDBhandle)
  
  # check if any bundles exist
  if(nrow(bndls) == 0){
    return()
  }
  
  progress = 0
  if(verbose){
    bundleCount=nrow(bndls)
    cat("INFO: Rewriting", bundleCount, "_annot.json files to file system...\n")
    pb=utils::txtProgressBar(min=0,max=bundleCount,style=3)
    utils::setTxtProgressBar(pb,progress)
  }
  
  for(i in 1:nrow(bndls)){
    bndl = bndls[i,]
    bundleAnnotDFs = load_bundleAnnotDFsDBI(emuDBhandle, bndl$session, bndl$name)
    annotJSONchar = bundleAnnotDFsToAnnotJSONchar(emuDBhandle, bundleAnnotDFs)
    
    # construct path to annotJSON
    annotFilePath = file.path(emuDBhandle$basePath, paste0(bndl$session, session.suffix), 
                              paste0(bndl$name, bundle.dir.suffix), 
                              paste0(bndl$name, bundle.annotation.suffix, '.json'))
    
    writeLines(annotJSONchar, annotFilePath)
    
    progress=progress+1L
    if(verbose){
      utils::setTxtProgressBar(pb,progress)
    }
  } 
}



#########################################################
# store / create / load functions

## Store EMU database to directory
## 
## @details 
## options is a list of key value pairs:
## rewriteSSFFTracks if TRUE rewrite SSF tracks instead of file copy to get rid of big endian encoded SSFF files (SPARC), default: FALSE
## ignoreMissingSSFFTrackFiles if TRUE missing SSFF track files are ignored, default: FALSE
## symbolicLinkSignalFiles if TRUE signal files are symbolic linked instead of copied. Implies: rewriteSSFFTracks=FALSE, Default: FALSE
## 
## @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
## @param targetDir target directory
## @param options list of options
## @param verbose show infos and progress bar
## @import stringr uuid jsonlite
## @keywords emuDB database Emu
## @seealso  \code{\link{load_emuDB}}
## @examples
## \dontrun{
## # Store database 'ae' to directory /homes/mylogin/EMUnew/
## 
##   store('ae',"/homes/mylogin/EmuStore/")
## 
## }
## 
##' @import stringr uuid jsonlite
store<-function(emuDBhandle, targetDir, options=NULL, verbose=TRUE){
  # default options
  # ignore missing SSFF track files
  # rewrite SSFF track files
  mergedOptions=list(ignoreMissingSSFFTrackFiles=TRUE,rewriteSSFFTracks=FALSE,symbolicLinkSignalFiles=FALSE)
  if(!is.null(options)){
    for(opt in names(options)){
      mergedOptions[[opt]]=options[[opt]]
    }
  }
  
  progress=0
  # check target dir
  if(file.exists(targetDir)){
    tdInfo=file.info(targetDir)
    if(!tdInfo[['isdir']]){
      stop(targetDir," exists and is not a directory.")
    }
  }else{
    # create target dir
    dir.create(targetDir)
  }
  
  # build db dir name
  dbDirName=paste0(emuDBhandle$dbName,emuDB.suffix)
  # create database dir in targetdir
  pp=file.path(targetDir,dbDirName)
  # check existence
  if(file.exists(pp)){
    stop(pp," already exists.")
  }
  
  dir.create(pp)
  
  # check if handle has basePath if not -> emuDB doesn't extist yet -> create new DBconfig
  if(is.null(emuDBhandle$basePath)){
    DBconfig = list(name = emuDBhandle$dbName, UUID=emuDBhandle$UUID, mediafileExtension = "wav", ssffTrackDefinitions=list(),levelDefinitions=list(),linkDefinitions=list())
  }else{
    DBconfig = load_DBconfig(emuDBhandle)
  }
  
  # set editable + showHierarchy
  DBconfig[['EMUwebAppConfig']][['activeButtons']]=list(saveBundle=TRUE,
                                                        showHierarchy=TRUE)
  
  # store db schema file
  store_DBconfig(emuDBhandle,DBconfig, basePath=pp)
  
  # create session dirs
  sessions = list_sessionsDBI(emuDBhandle)
  if(nrow(sessions) == 0){
    return()
  }
  sesDirPaths = file.path(pp, paste0(sessions$name, session.suffix))
  for(path in sesDirPaths){
    dir.create(path)
  }
  
  # create bundle dirs
  bndls = list_bundlesDBI(emuDBhandle)
  if(nrow(bndls) == 0){
    return()
  }
  bndlDirPaths = file.path(pp, paste0(sessions$name, session.suffix), paste0(bndls$name, bundle.dir.suffix))
  for(path in bndlDirPaths){
    dir.create(path)
  }
  
  # copy media files
  mediaFilePathsOld = file.path(emuDBhandle$basePath, paste0(sessions$name, session.suffix), paste0(bndls$name, bundle.dir.suffix), paste0(bndls$name, ".", DBconfig$mediafileExtension))
  mediaFilePathsNew = file.path(pp, paste0(sessions$name, session.suffix), paste0(bndls$name, bundle.dir.suffix), paste0(bndls$name, ".", DBconfig$mediafileExtension))
  file.copy(mediaFilePathsOld, mediaFilePathsNew)
  
  # rewrite annotations (or should these just be a copied as well?)
  emuDBhandle$basePath = pp
  rewrite_allAnnots(emuDBhandle, verbose = verbose)
  
  # copy SSFF files
  ssffDefs = list_ssffTrackDefinitions(emuDBhandle)
  if(!is.null(ssffDefs)){
    for(i in 1:nrow(ssffDefs)){
      ssffDef = ssffDefs[1,]
      ssffFilePathsOld = file.path(emuDBhandle$basePath, paste0(sessions$name, session.suffix), paste0(bndls$name, bundle.dir.suffix), paste0(bndls$name, ".", ssffDef$fileExtension))
      ssffFilePathsNew = file.path(pp, paste0(sessions$name, session.suffix), paste0(bndls$name, bundle.dir.suffix), paste0(bndls$name, ".", ssffDef$fileExtension))
      file.copy(ssffFilePathsOld, ssffFilePathsNew)
    }
  }
  
}


##' @title Create empty emuDB
##' @description Creates an empty emuDB in the target directory specified
##' @details Creates a new directory [name]_emuDB in targetDir. By default the emuDB is created in the R session, 
##' written to the filesystem and then purged from the R session.
##' @param name of new emuDB
##' @param targetDir target directory to store the emuDB to
##' @param mediaFileExtension defines mediaFileExtention (NOTE: currently only 
##' 'wav' (the default) is supported by all components of EMU)
##' @param store store new created emuDB to file system
##' @param verbose display infos & show progress bar
##' @export
##' @examples 
##' \dontrun{
##' # create empty emuDB in folder provided by tempdir()
##' create_emuDB(name = "myNewEmuDB", 
##'              targetDir = tempdir())
##' }
create_emuDB<-function(name, targetDir, mediaFileExtension='wav', 
                       store=TRUE, verbose=TRUE){
  
  dbDirName=paste0(name,emuDB.suffix)
  dbHandle = emuDBhandle(dbName = name , basePath=NULL, uuid::UUIDgenerate(), ":memory:")
  if(store){
    store(dbHandle, targetDir=targetDir, verbose = verbose)
  }
  
  return(invisible())
}

##' Load emuDB
##' 
##' @description Function loads emuDB into its cached representation and makes it accessible from within the 
##' current R session by returning a emuDBhandle object
##' @details In order to access an emuDB from R it is necessary to load the annotation and configuration 
##' files to an emuR internal database format. The function expects a emuDB file structure in directory 
##' \code{databaseDir}. The emuDB configuration file is loaded first. On success the function iterates 
##' through session and bundle directories and loads found annotation files. The parameter \code{inMemoryCache} 
##' determines where the internal database is stored: If \code{FALSE} a databse cache file in \code{databaseDir} 
##' is used. When the database is loaded for the first time the function will create a new cache file and store 
##' the data to it. On subsequent loading of the same database the cache is only updated if files have changed, 
##' therefore the loading is then much faster. For this to work the user needs write permissions to 
##' \code{databaseDir} and the cache file. The database is loaded into a volatile in-memory database if 
##' \code{inMemoryCache} is set to \code{TRUE}.
##' @param databaseDir directory of the emuDB
##' @param inMemoryCache cache the loaded DB in memory
##' @param connection pass in DBI connection to SQL database if you want to override the default which is to 
##' use an SQLite database either in memory (\code{inMemoryCache = TRUE}) or in the emuDB folder. This is intended
##' for expert use only!
##' @param verbose be verbose
##' @param ... additional parameters
##' @return name of emuDB
##' @import jsonlite DBI
##' @export
##' @keywords emuDB database DBconfig
##' @examples
##' \dontrun{
##' ## Load database ae in directory /homes/mylogin/EMUnew/ae 
##' ## assuming an existing emuDB structure in this directory
##' 
##' ae = load_emuDB("/homes/mylogin/EMU/ae")
##' 
##' ## Load database ae from demo data
##' 
##' # create demo data in temporary directory
##' create_emuRdemoData(dir = tempdir())
##' # build base path to demo emuDB
##' demoDatabaseDir = file.path(tempdir(), "emuR_demoData", "ae_emuDB")
##' 
##' # load demo emuDB
##' ae = load_emuDB(demoDatabaseDir)
##' 
##' }
load_emuDB <- function(databaseDir, inMemoryCache = FALSE, connection = NULL, verbose=TRUE, ...){
  progress = 0
  # check database dir
  if(!file.exists(databaseDir)){
    stop("Database dir ",databaseDir," does not exist!")
  }
  dbDirInfo=file.info(databaseDir)
  if(!dbDirInfo[['isdir']]){
    stop(databaseDir," exists, but is not a directory.")
  }
  
  # load db schema file
  dbCfgPattern=paste0('.*',database.schema.suffix,'$')
  dbCfgFiles=list.files(path=databaseDir,dbCfgPattern)
  dbCfgFileCount=length(dbCfgFiles)
  if(dbCfgFileCount==0){
    stop("Could not find global DB config JSON file (regex pattern: ",dbCfgPattern,") in ",databaseDir)
  }
  if(dbCfgFileCount>1){
    stop("Found multiple global DB config JSON files (regex pattern: ",dbCfgPattern,") in ",databaseDir)
  }
  
  dbCfgPath=file.path(databaseDir,dbCfgFiles[[1]])
  if(!file.exists(dbCfgPath)){
    stop("Could not find database info file: ",dbCfgPath,"\n")
  }
  # extract ... (ellipsis) parameters
  dots = list(...)
  if("update_cache" %in% names(dots)){
    updateCache = dots$update_cache
  }else{
    updateCache = T
  }
  
  # load DBconfig
  DBconfig = jsonlite::fromJSON(dbCfgPath, simplifyVector=FALSE)
  # normalize base path
  basePath = normalizePath(databaseDir)
  
  # shorthand vars
  dbName = DBconfig$name
  dbUUID = DBconfig$UUID
  
  # create dbHandle
  if(inMemoryCache){
    dbHandle = emuDBhandle(dbName, basePath, dbUUID, connectionPath = ":memory:")
  }else{
    cachePath = file.path(normalizePath(databaseDir), paste0(dbName, database.cache.suffix))
    if(is.null(connection)){
      dbHandle = emuDBhandle(dbName, basePath, dbUUID, cachePath)
    }else{
      dbHandle = emuDBhandle(dbName, basePath, dbUUID, "", connection = connection)
    }
  }
  # check if cache exist -> update cache if true
  dbsDf = get_emuDbDBI(dbHandle)
  if(nrow(dbsDf)>0){
    if(updateCache){
      update_cache(dbHandle, verbose = verbose)
    }
    return(dbHandle)
  }
  
  # write to DBI emuDB table
  add_emuDbDBI(dbHandle)
  
  # list sessions & bundles
  sessions = list_sessions(dbHandle)
  bundles = list_bundles(dbHandle)
  # add column to sessions to track if already stored
  if(nrow(sessions) != 0){
    sessions$stored = F
    
    # calculate bundle count
    bundleCount = nrow(bundles)
    # create progress bar
    pMax = bundleCount
    if(pMax == 0){
      pMax = 1
    }
    if(verbose){ 
      cat(paste0("INFO: Loading EMU database from ", databaseDir, "... (", bundleCount , " bundles found)\n"))
      pb=utils::txtProgressBar(min = 0L, max = pMax, style = 3)
      utils::setTxtProgressBar(pb, progress)
    }
    # bundles
    for(bndlIdx in 1:nrow(bundles)){
      bndl = bundles[bndlIdx,]
      # check if session has to be added to DBI
      if(!(sessions$stored[sessions$name == bndl$session])){
        add_sessionDBI(dbHandle, bndl$session)
        sessions$stored[sessions$name == bndl$session] = TRUE
      }
      
      # construct path to annotJSON
      annotFilePath = normalizePath(file.path(dbHandle$basePath, paste0(bndl$session, session.suffix), 
                                              paste0(bndl$name, bundle.dir.suffix), 
                                              paste0(bndl$name, bundle.annotation.suffix, '.json')))
      
      # calculate MD5 sum of bundle annotJSON
      newMD5annotJSON = tools::md5sum(annotFilePath)
      names(newMD5annotJSON) = NULL
      
      # read annotJSON as charac 
      annotJSONchar = readChar(annotFilePath, file.info(annotFilePath)$size)
      
      # convert to bundleAnnotDFs
      bundleAnnotDFs = annotJSONcharToBundleAnnotDFs(annotJSONchar)
      
      # add to bundle table
      add_bundleDBI(dbHandle, bndl$session, bndl$name, bundleAnnotDFs$annotates, bundleAnnotDFs$sampleRate, newMD5annotJSON)
      # add to items, links, labels tables
      store_bundleAnnotDFsDBI(dbHandle, bundleAnnotDFs, bndl$session, bndl$name)
      
      # increase progress bar  
      progress=progress+1L
      if(verbose){
        utils::setTxtProgressBar(pb,progress)
      }
      
    }
    
    # build redundat links and calc positions
    if(verbose){ 
      cat("\nbuilding redundant links and position of links... (this may take a while)\n")
    }else{
      cat("\n")
    }
    
    build_allRedundantLinks(dbHandle)
    calculate_postionsOfLinks(dbHandle)
  }
  
  return(dbHandle)
  
}

#######################
# FOR DEVELOPMENT
# library('testthat')
# test_file('tests/testthat/test_aaa_initData.R')
# test_file('tests/testthat/test_emuR-database.R')
# test_file('tests/testthat/test_duplicate.loaded.emuDB.R')
# test_file('tests/testthat/test_database.caching.R')