
convert_queryResultToEmusegs<-function(emuDBhandle, timeRefSegmentLevel=NULL, filteredTablesSuffix){
  emuRsegs = convert_queryResultToEmuRsegs(emuDBhandle, timeRefSegmentLevel, filteredTablesSuffix)
  emusegs = as.emusegs(emuRsegs)
  return(emusegs)
}

##################################
#
fconvert_queryResultToEmuRsegs <- function(emuDBhandle, timeRefSegmentLevel=NULL, filteredTablesSuffix){
  itemsTableName = paste0("items", filteredTablesSuffix)
  labelsTableName = paste0("labels", filteredTablesSuffix)
  linksExtTableName = paste0("links_ext", filteredTablesSuffix)
  
  projectionItemsN = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT COUNT(*) AS n FROM interm_res_proj_items_tmp_root"))$n
  if(projectionItemsN > 0){ 
    # use projection items
    itsTableName = "interm_res_proj_items_tmp_root"
    seqStartIdColName = "p_seq_start_id"
    seqEndIdColName = "p_seq_end_id"
    seqLenColName = "p_seq_len"
    levelColName = "p_level"
  }else{
    # use "normal" items
    itsTableName = "interm_res_items_tmp_root"
    seqStartIdColName = "seq_start_id"
    seqEndIdColName = "seq_end_id"
    seqLenColName = "seq_len"
    levelColName = "level"
  }
  
  
  dbConfig = load_DBconfig(emuDBhandle)
  levelNamesWithTime = unlist(lapply(dbConfig$levelDefinitions, function(x){if(x$type=="SEGMENT" || x$type=="EVENT"){return(x$name)} }))
  
  resultLevel = DBI::dbGetQuery(emuDBhandle$connection, "SELECT * FROM interm_res_meta_infos_tmp_root")$result_level
  if(resultLevel %in% levelNamesWithTime){
    
    seglist = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT labels.label AS labels, ",
                                                             "CASE items_seq_start.type ",
                                                             " WHEN 'SEGMENT' THEN items_seq_start.sample_start / items_seq_start.sample_rate * 1000.0 ",
                                                             " WHEN 'EVENT' THEN 'NOT IMPLEMENTED YET' ",
                                                             " ELSE 'SIC!! Something went wrong' ",
                                                             "END AS start, ",
                                                             "CASE items_seq_start.type ",
                                                             " WHEN 'SEGMENT' THEN (items_seq_end.sample_start + items_seq_end.sample_dur) / items_seq_end.sample_rate * 1000.0 ",
                                                             " WHEN 'EVENT' THEN 'NOT IMPLEMENTED YET' ",
                                                             " ELSE 'SIC!! Something went wrong' ",
                                                             "END AS end, ",
                                                             "interm_res_items_tmp_root.session || ':' || interm_res_items_tmp_root.bundle AS utts, ",
                                                             "interm_res_items_tmp_root.db_uuid, interm_res_items_tmp_root.session, interm_res_items_tmp_root.bundle, interm_res_items_tmp_root.seq_start_id AS startItemID, interm_res_items_tmp_root.seq_end_id AS endItemID, ",
                                                             "interm_res_items_tmp_root.level AS level, items_seq_start.type AS type, ",
                                                             "items_seq_start.sample_start AS sampleStart, (items_seq_end.sample_start + items_seq_end.sample_dur) AS sampleEnd, items_seq_start.sample_rate AS sampleRate ",
                                                             "FROM interm_res_items_tmp_root, items AS items_seq_start, items AS items_seq_end, labels ",
                                                             "WHERE interm_res_items_tmp_root.db_uuid = items_seq_start.db_uuid AND interm_res_items_tmp_root.session = items_seq_start.session AND interm_res_items_tmp_root.bundle = items_seq_start.bundle AND interm_res_items_tmp_root.seq_start_id = items_seq_start.item_id ",
                                                             "AND interm_res_items_tmp_root.db_uuid = items_seq_end.db_uuid AND interm_res_items_tmp_root.session = items_seq_end.session AND interm_res_items_tmp_root.bundle = items_seq_end.bundle AND interm_res_items_tmp_root.seq_end_id = items_seq_end.item_id ",
                                                             "AND interm_res_items_tmp_root.db_uuid = labels.db_uuid AND interm_res_items_tmp_root.session = labels.session AND interm_res_items_tmp_root.bundle = labels.bundle AND interm_res_items_tmp_root.seq_end_id = labels.item_id"))
    
  }else{
    
    # for(lnwt in levelNamesWithTime){
    lnwt = levelNamesWithTime[1]
    connectHierPaths = get_hierPathsConnectingLevels(emuDBhandle, lnwt, resultLevel)
    
    tln = connectHierPaths[[1]][length(connectHierPaths[[1]])]
    # insert all time items into new table
    timeItemsTableSuffix = "time_level_items"
    create_intermResTmpQueryTablesDBI(emuDBhandle, suffix = timeItemsTableSuffix)
    DBI::dbGetQuery(emuDBhandle$connection, paste0("INSERT INTO interm_res_items_tmp_", timeItemsTableSuffix, " ",
                                                   "SELECT db_uuid, session, bundle, item_id AS seq_start_id, item_id AS seq_end_id, 1 AS seq_len, level  FROM ", itemsTableName, " ",
                                                   "WHERE db_uuid ='", emuDBhandle$UUID, "' AND level = '", tln, "'"))
    
    query_databaseHier(emuDBhandle, tln, resultLevel, timeItemsTableSuffix, "root", filteredTablesSuffix) # result written to lr_exp_res_tmp table
    
    # seglist=DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT \
    #                      labels,
    #                                                        CASE type WHEN 'EVENT' THEN \
    #                                                        CAST (sample_start AS REAL)/ CAST( sample_rate AS REAL) * 1000.0 \
    #                                                        ELSE \
    #                                                        (CAST (sample_start AS REAL) + 0.5 ) / CAST( sample_rate AS REAL) * 1000.0 \
    #                                                        END AS start, \
    #                                                        CASE type WHEN 'EVENT' THEN \
    #                                                        0.0
    #                                                        ELSE \
    #                                                        (CAST (sample_end AS REAL) + 1.5 ) / CAST( sample_rate AS REAL) * 1000.0 \
    #                                                        END AS end, \
    #                                                        session || ':' || bundle AS utts, \
    #                                                        db_uuid,session,bundle, start_item_id AS startItemID, end_item_id AS endItemID,level,type,sample_start AS sampleStart,sample_end AS sampleEnd,sample_rate AS sampleRate\
    #                                                        FROM (", queryStr, ")"))
    seglist = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT 'XXX' AS labels, ",
                                                             "CASE type ",
                                                             " WHEN 'SEGMENT' THEN min(items.sample_start) / items.sample_rate * 1000.0 ",
                                                             " ELSE 'SIC SIC SIC' ",
                                                             "END AS start, ",
                                                             "CASE type ",
                                                             " WHEN 'SEGMENT' THEN max(items.sample_start + items.sample_dur) / items.sample_rate * 1000.0 ",
                                                             " ELSE 'SIC SIC SIC' ",
                                                             "END AS end, ",
                                                             "items.session || ':' || items.bundle AS utts, ",
                                                             "items.db_uuid, items.session, items.bundle, lr_exp_res_tmp.r_seq_start_id AS startItemID, lr_exp_res_tmp.r_seq_end_id AS endItemID, ",
                                                             "'todo' AS level, 'TODO' AS type, ",
                                                             "min(items.sample_start) AS sampleStart, max(items.sample_start + items.sample_dur) AS sampleEnd, items.sample_rate AS sampleRate ",
                                                             "FROM lr_exp_res_tmp, items", filteredTablesSuffix, " ",
                                                             "WHERE lr_exp_res_tmp.db_uuid = items.db_uuid AND lr_exp_res_tmp.session = items.session AND lr_exp_res_tmp.bundle = items.bundle AND lr_exp_res_tmp.l_seq_start_id = items.item_id ",
                                                             "GROUP BY r_seq_start_id"))
  }
  # set emusegs type attribute, default 'segment'
  slType='segment'
  if(nrow(seglist)>0){
    # set to event only if all rows are of type EVENT
    if(all(seglist$type == "EVENT")){
      slType='event'
    }
  }
  queryStr = DBI::dbGetQuery(emuDBhandle$connection, "SELECT query_str FROM interm_res_meta_infos_tmp_root")$query_str
  segmentList=make.emuRsegs(dbName = emuDBhandle$dbName, seglist = seglist, query = queryStr, type = slType)
  return(segmentList)
  
}

#################################
#
convert_queryResultToEmuRsegs <- function(emuDBhandle, timeRefSegmentLevel=NULL, filteredTablesSuffix){
  itemsTableName = paste0("items", filteredTablesSuffix)
  labelsTableName = paste0("labels", filteredTablesSuffix)
  linksExtTableName = paste0("links_ext", filteredTablesSuffix)
  
  bundles=c()
  labels=c()
  start=c()
  end=c()
  slType=NULL
  projectionItemsN = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT COUNT(*) AS n FROM interm_res_proj_items_tmp_root"))$n
  if(projectionItemsN > 0){ 
    # use projection items
    itsTableName = "interm_res_proj_items_tmp_root"
    seqStartIdColName = "p_seq_start_id"
    seqEndIdColName = "p_seq_end_id"
    seqLenColName = "p_seq_len"
    levelColName = "p_level"
  }else{
    # use "normal" items
    itsTableName = "interm_res_items_tmp_root"
    seqStartIdColName = "seq_start_id"
    seqEndIdColName = "seq_end_id"
    seqLenColName = "seq_len"
    levelColName = "level"
  }
  
  # get distinct result levels
  distinctLevels = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT DISTINCT ", levelColName, " FROM ", itsTableName))
  for(attrNm in distinctLevels[,levelColName]){
    
    lvlNm = get_levelNameForAttributeName(emuDBhandle, attributeName = attrNm)
    ld = get_levelDefinition(emuDBhandle, name = lvlNm)
    
    if(ld['type']=='ITEM'){
      segLvlNms = find_segmentLevels(emuDBhandle, attrNm)
      if(!is.null(timeRefSegmentLevel)){
        if(!(timeRefSegmentLevel %in% segLvlNms)){
          stop("Cannot resolve time information for result level '",attrNm,"' using segment time reference level '",timeRefSegmentLevel,"'\nPlease set one of these levels for timeRefSegmentLevel parameter: ",paste(segLvlNms,collapse=', '),".")
        }
      }else{
        segLvlsCnt=length(segLvlNms)
        if(segLvlsCnt>1){
          
          stop("Segment time information derivation for level '",attrNm,"' is ambiguous:\nThe level is linked to multiple segment levels: ",paste(segLvlNms,collapse=', '),"\nPlease select one of these levels using the 'timeRefSegmentLevel' query parameter.")
        }
      }
    }
    
  }
  
  itCount = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT COUNT(*) AS n FROM ", itsTableName))$n
  if(itCount > 0){
    maxSeqLen = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT max(", seqLenColName, ") AS max_seq_len FROM ", itsTableName))$max_seq_len
  }else{
    maxSeqLen=1L
  }
  
  # query seglist data except labels
  # for this data the information in start end item of the sequence is sufficient
  # it takes only the start and end items of the query result into account
  # the CASE WHEN THEN ELSE END terms are necessary to get the start and end samples of sequences which are not segment levels and therefore have no time information
  linksExtFilteredCount = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT COUNT(*) AS n FROM ", linksExtTableName))$n
  hasLinks = (linksExtFilteredCount > 0)
  
  # select columns: id,session,bundle,start_item_id,end_item_id ,type ...
  selectStr=paste0("SELECT s.db_uuid ,s.session,s.bundle,s.item_id AS start_item_id ,e.item_id AS end_item_id,r.", levelColName, ",s.type, ")
  
  # find sequence start position
  # use sample start of sequence start item for type SEGMENT and samplePoint for type EVENT 
  selectStr=paste0(selectStr,"CASE s.type \
                   WHEN 'SEGMENT' THEN s.sample_start \
                   WHEN 'EVENT' THEN s.sample_point ")
  # calculate start sample for type ITEM
  if(hasLinks){
    # items of type ITEM have no (sample) time information
    # therefore we search for linked SEGMENT items and take their start sample position
    selectStr=paste0(selectStr," ELSE (SELECT i.sample_start FROM ", itemsTableName, " i WHERE i.db_uuid=s.db_uuid AND i.session=s.session AND i.bundle=s.bundle AND i.type='SEGMENT' AND ")
    if(!is.null(timeRefSegmentLevel)){
      selectStr=paste0(selectStr," i.level='",timeRefSegmentLevel,"' AND ")
    }
    selectStr=paste0(selectStr," EXISTS (SELECT * FROM ", linksExtTableName, " l WHERE s.db_uuid=s.db_uuid AND s.session=l.session AND s.bundle=l.bundle AND s.item_id=l.from_id AND i.db_uuid=l.db_uuid AND i.session=l.session AND i.bundle=l.bundle AND i.item_id=l.to_id AND l.to_seq_idx=0)) ")
  }else{
    # TODO
    # No sample start information. (throw error ?)
  }
  selectStr=paste0(selectStr," END AS sample_start, ")
  
  # find sequence end position
  # use sample start plus sample duration of sequence end item for type SEGMENT and zero for type EVENT
  # TODO is zero correct here ?? should it be e.sample_point instead ??
  selectStr=paste0(selectStr," CASE s.type \
                   WHEN 'SEGMENT' THEN (e.sample_start+e.sample_dur) \
                   WHEN 'EVENT' THEN NULL ")  
  if(hasLinks){
    # items of type ITEM have no (sample) time information
    # therefore we search for linked SEGMENT items and take their end sample position
    selectStr=paste0(selectStr," ELSE (SELECT i.sample_start+i.sample_dur FROM ", itemsTableName, " i WHERE i.db_uuid=e.db_uuid AND i.session=e.session AND i.bundle=e.bundle AND i.type='SEGMENT' AND ")
    if(!is.null(timeRefSegmentLevel)){
      selectStr=paste0(selectStr," i.level='",timeRefSegmentLevel,"' AND ")
    }
    selectStr=paste0(selectStr," EXISTS (SELECT * FROM ", linksExtTableName, " l WHERE e.db_uuid=l.db_uuid AND e.session=l.session AND e.bundle=l.bundle AND e.item_id=l.from_id AND i.db_uuid=l.db_uuid AND i.session=l.session AND i.bundle=l.bundle AND i.item_id=l.to_id AND l.to_seq_idx+1=l.to_seq_len)) ")
    
  }
  selectStr=paste0(selectStr,"END AS sample_end, ")
  
  # find samplerate
  # use sample rate of sequence start item for type SEGMENT and EVENT 
  selectStr=paste0(selectStr,"CASE s.type \
                   WHEN 'SEGMENT' THEN s.sample_rate \
                   WHEN 'EVENT' THEN s.sample_rate ")
  if(hasLinks){
    # items of type ITEM have no sample rate information
    # therefore we search for linked SEGMENT items and take their start sample position
    # TODO Can we use EVENT items as well ?
    selectStr=paste0(selectStr," ELSE (SELECT i.sample_rate FROM ", itemsTableName, " i WHERE i.db_uuid=s.db_uuid AND i.session=s.session AND i.bundle=s.bundle AND i.type='SEGMENT' AND ")
    if(!is.null(timeRefSegmentLevel)){
      selectStr=paste0(selectStr," i.level='",timeRefSegmentLevel,"' AND ")
    }
    selectStr=paste0(selectStr," EXISTS (SELECT * FROM ", linksExtTableName, " l WHERE s.item_id=l.from_id AND i.item_id=l.to_id AND i.db_uuid=l.db_uuid AND i.session=l.session AND i.bundle=l.bundle AND l.to_seq_idx=0)) ")
    
  }else{
    # TODO no samplerate , error ??
  }
  
  selectStr=paste0(selectStr," END AS sample_rate, ")
  
  # from clause
  fromStr=paste0("FROM ", itemsTableName, " s,", itemsTableName, " e,", itsTableName, " r, ")
  
  # where clause: make sure start and end are in same emuDB, session and bundle, select start and end id
  whereStr=paste0("WHERE e.db_uuid=s.db_uuid AND e.session=s.session AND e.bundle=s.bundle AND r.db_uuid=s.db_uuid AND r.session=s.session AND r.bundle=s.bundle AND s.item_id=", seqStartIdColName, " AND e.item_id=", seqEndIdColName ," AND e.level=s.level AND ")
  
  # order
  orderStr="ORDER BY s.db_uuid,s.session,s.bundle,start_item_id,end_item_id"
  
  # append terms depending on maximum sequence length
  # build query for label sequence string
  # TODO assume equal seq len for each sequence for now!!
  # this is not the case for queries like emu.query("andosl","*","[#[Phonetic=t -> Phonetic=S] -> #Phonetic=I]")
  # which are not allowed by BNF and in emuR but in fact they are working with Emu 2.3 and emuR requery may produce such results
  # result would be a mix with t->S and I items (sequence lengths 2 and 1)
  for(seqIdx in 1:maxSeqLen){
    selectStr=paste0(selectStr, "(SELECT l.label FROM ", labelsTableName, " l WHERE l.db_uuid=i", seqIdx, ".db_uuid AND l.session=i", seqIdx, ".session AND l.bundle=i", seqIdx, ".bundle AND l.item_id=i", seqIdx, ".item_id AND l.name=r.", levelColName, ")")
    
    fromStr=paste0(fromStr, itemsTableName, " i",seqIdx)
    offset=seqIdx-1
    whereStr=paste0(whereStr, "i", seqIdx, ".db_uuid=s.db_uuid AND i", seqIdx, ".session=s.session AND i", seqIdx, ".bundle=s.bundle AND i", seqIdx ,".level=s.level AND i", seqIdx, ".seq_idx=s.seq_idx+", offset)
    if(seqIdx<maxSeqLen){
      selectStr=paste0(selectStr," || '->' || ")
      fromStr=paste0(fromStr,',')
      whereStr=paste0(whereStr,' AND ')
    }
  }
  selectStr=paste0(selectStr," AS labels ")
  
  queryStr=paste(selectStr,fromStr,whereStr,orderStr,sep = ' ')
  
  # convert samples to milliseconds SQL string:
  queryStrInclConvert = paste0("SELECT labels,
                       CASE type WHEN 'EVENT' THEN \
                        CAST (sample_start AS REAL)/ CAST( sample_rate AS REAL) * 1000.0 \
                       ELSE \
                        (CAST (sample_start AS REAL) + 0.5 ) / CAST( sample_rate AS REAL) * 1000.0 \
                       END AS start, \
                       CASE type WHEN 'EVENT' THEN \
                         0.0
                       ELSE \
                        (CAST (sample_end AS REAL) + 1.5 ) / CAST( sample_rate AS REAL) * 1000.0 \
                       END AS end, \
                       session || ':' || bundle AS utts, \
                       db_uuid,session,bundle, start_item_id  AS startItemID, end_item_id AS endItemID,", levelColName, " AS level,type, sample_start AS sampleStart,sample_end AS sampleEnd,sample_rate AS sampleRate \
                      FROM (", queryStr, ") AS qs_res ORDER BY db_uuid,session,bundle,sampleStart")
  
  seglist = DBI::dbGetQuery(emuDBhandle$connection, queryStrInclConvert)
  
  # set emusegs type attribute, default 'segment'
  slType='segment'
  if(nrow(seglist)>0){
    # set to event only if all rows are of type EVENT
    if(all(seglist$type == "EVENT")){
      slType='event'
    }
  }
  queryStr = DBI::dbGetQuery(emuDBhandle$connection, "SELECT query_str FROM interm_res_meta_infos_tmp_root")$query_str
  segmentList=make.emuRsegs(dbName = emuDBhandle$dbName, seglist = seglist, query = queryStr, type = slType)
  return(segmentList)
}


##################################
##################################
##################################
# convert to emuRsegs segemnt list, variable sequenec length of input allowed
convert_queryResultToVariableEmuRsegs <- function(emuDBhandle, timeRefSegmentLevel=NULL){
  its=NULL
  
  bundles=c()
  labels=c()
  start=c()
  end=c()
  slType=NULL
  
  itsTableName = "interm_res_items_tmp_root"
  
  # get distinct result levels
  distinctLevels=DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT DISTINCT level FROM ", itsTableName))
  for(attrNm in distinctLevels[,'level']){
    lvlNm=get_levelNameForAttributeName(emuDBhandle,attributeName = attrNm)
    ld=get_levelDefinition(emuDBhandle, name = lvlNm)
    
    if(ld['type']=='ITEM'){
      segLvlNms=find_segmentLevels(emuDBhandle,attrNm)
      if(!is.null(timeRefSegmentLevel)){
        if(!(timeRefSegmentLevel %in% segLvlNms)){
          stop("Cannot resolve time information for result level '",attrNm,"' using segment time reference level '",timeRefSegmentLevel,"'\nPlease set one of these levels for timeRefSegmentLevel parameter: ",paste(segLvlNms,collapse=', '),".")
        }
      }else{
        segLvlsCnt=length(segLvlNms)
        if(segLvlsCnt>1){
          
          stop("Segment time information derivation for level '",attrNm,"' is ambiguous:\nThe level is linked to multiple segment levels: ",paste(segLvlNms,collapse=', '),"\nPlease select one of these levels using the 'timeRefSegmentLevel' query parameter.")
        }
      }
    }
    
  }
  
  
  # itCount=nrow(its)
  itCount = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT COUNT(*) AS N FROM ", itsTableName))$N
  if(itCount==0){
    its=data.frame(db_uuid=character(0),session=character(0),bundle=character(0),seq_start_id=integer(0),seq_end_id=integer(0),seq_len=integer(0),level=character(0),stringsAsFactors = FALSE)
  }
  
  if(itCount>0){
    maxSeqLenDf=DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT max(seq_len) AS maxSeqLen FROM ", itsTableName))
    maxSeqLen=maxSeqLenDf[1,'maxSeqLen']
    
    # for string conacatenation: we need all occuring seq lengths 
    # distinct sequence lengths
    # distinctSeqLens=sqldf(c(resIdxSql,"SELECT DISTINCT seq_len FROM its"))
    distinctSeqLens = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT DISTINCT seq_len FROM ", itsTableName))
  }else{
    maxSeqLen=1L
  }
  
  
  # query seglist data except labels
  # for this data the information in start end item of the sequence is sufficient
  # it takes only the start  and end items of the query result in account
  # the CASE WHEN THEN ELSE END terms are necessary to get the start and end samples of sequences which are not segment levels and therefore have no time information  
  hasLinks = DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT COUNT(*) AS n FROM links_ext"))$n > 0
  
  
  # select columns: id,session,bundle,start_item_id,end_item_id ,type ...
  selectStr=paste0("SELECT s.db_uuid ,s.session,s.bundle,s.item_id AS start_item_id ,e.item_id AS end_item_id, r.level, s.type, ")
  
  # find sequence start position
  # use sample start of sequence start item for type SEGMENT and sample_point for type EVENT 
  selectStr=paste0(selectStr,"CASE s.type \
                   WHEN 'SEGMENT' THEN s.sample_start \
                   WHEN 'EVENT' THEN s.sample_point ")
  # calculate start sample for type ITEM
  if(hasLinks){
    # items of type ITEM have no (sample) time information
    # therefore we search for linked SEGMENT items and take their start sample position
    selectStr=paste0(selectStr," ELSE (SELECT i.sample_start FROM items i WHERE i.db_uuid=s.db_uuid AND i.session=s.session AND i.bundle=s.bundle AND i.type='SEGMENT' AND ")
    if(!is.null(timeRefSegmentLevel)){
      selectStr=paste0(selectStr," i.level='",timeRefSegmentLevel,"' AND ")
    }
    selectStr=paste0(selectStr," EXISTS (SELECT * FROM links_ext l WHERE s.db_uuid=s.db_uuid AND s.session=l.session AND s.bundle=l.bundle AND s.item_id=l.from_id AND i.db_uuid=l.db_uuid AND i.session=l.session AND i.bundle=l.bundle AND i.item_id=l.to_id AND l.to_seq_idx=0)) ")
  }else{
    # TODO
    # No sample start information. (throw error ?)
  }
  selectStr=paste0(selectStr," END AS sample_start, ")
  
  # find sequence end position
  # use sample start plus sample duration of sequence end item for type SEGMENT and zero for type EVENT
  # TODO is zero correct here ?? should it be e.samplePoint instead ??
  selectStr=paste0(selectStr," CASE s.type \
                   WHEN 'SEGMENT' THEN (e.sample_start+e.sample_dur) \
                   WHEN 'EVENT' THEN NULL ")  
  if(hasLinks){
    # items of type ITEM have no (sample) time information
    # therefore we search for linked SEGMENT items and take their end sample position
    
    selectStr=paste0(selectStr," ELSE (SELECT i.sample_start+i.sample_dur FROM items i WHERE i.db_uuid=e.db_uuid AND i.session=e.session AND i.bundle=e.bundle AND i.type='SEGMENT' AND ")
    if(!is.null(timeRefSegmentLevel)){
      selectStr=paste0(selectStr," i.level='",timeRefSegmentLevel,"' AND ")
    }
    selectStr=paste0(selectStr," EXISTS (SELECT * FROM links_ext l WHERE e.db_uuid=l.db_uuid AND e.session=l.session AND e.bundle=l.bundle AND e.item_id=l.from_id AND i.db_uuid=l.db_uuid AND i.session=l.session AND i.bundle=l.bundle AND i.item_id=l.to_id AND l.to_seq_idx+1=l.to_seq_len)) ")
    
  }
  selectStr=paste0(selectStr,"END AS sample_end, ")
  
  # find samplerate
  # use sample rate of sequence start item for type SEGMENT and EVENT 
  selectStr=paste0(selectStr,"CASE s.type \
                   WHEN 'SEGMENT' THEN s.sample_rate \
                   WHEN 'EVENT' THEN s.sample_rate ")
  if(hasLinks){
    # items of type ITEM have no sample rate information
    # therefore we search for linked SEGMENT items and take their start sample position
    # TODO Can we use EVENT items as well ?
    selectStr=paste0(selectStr," ELSE (SELECT i.sample_rate FROM items i WHERE i.db_uuid=s.db_uuid AND i.session=s.session AND i.bundle=s.bundle AND i.type='SEGMENT' AND ")
    if(!is.null(timeRefSegmentLevel)){
      selectStr=paste0(selectStr," i.level='",timeRefSegmentLevel,"' AND ")
    }
    selectStr=paste0(selectStr," EXISTS (SELECT * FROM links_ext l WHERE s.item_id=l.from_id AND i.item_id=l.to_id AND i.db_uuid=l.db_uuid AND i.session=l.session AND i.bundle=l.bundle AND l.to_seq_idx=0)) ")
    
  }else{
    # TODO no samplerate , error ??
  }
  
  selectStr=paste0(selectStr," END AS sample_rate, ")
  
  # from clause
  fromStr=paste0("FROM items s,items e,", itsTableName, " r ")
  
  # where clause: make sure start and end are in same emuDB, session and bundle, select start and end id
  whereStr=paste0("WHERE e.db_uuid=s.db_uuid AND e.session=s.session AND e.bundle=s.bundle AND r.db_uuid=s.db_uuid AND r.session=s.session AND r.bundle=s.bundle AND s.item_id=r.seq_start_id AND e.item_id=r.seq_end_id AND e.level=s.level ")
  
  # order
  orderStr="ORDER BY s.db_uuid,s.session,s.bundle,start_item_id,end_item_id"
  
  if(itCount>0){
    selectStr=paste0(selectStr," CASE r.seq_len ")
    distinctSeqLensLength=nrow(distinctSeqLens)
    # for each seq len, which occurs in the input segment list
    for(si in 1:distinctSeqLensLength){
      seqLen=distinctSeqLens[si,1]
      # append terms depending on maximum sequence length
      # build query for label sequence string
      # emuR hierachical requery produces results with differnet seq lengths per row
      selectStr=paste0(selectStr,"WHEN ",seqLen," THEN ")
      for(seqIdx in 1:seqLen){
        selectStr=paste0(selectStr,'(SELECT l.label FROM labels l,items i WHERE l.db_uuid=i.db_uuid AND l.session=i.session AND l.bundle=i.bundle AND l.item_id=i.item_id AND l.name=r.level AND ')
        
        offset=seqIdx-1
        selectStr=paste0(selectStr,'i.db_uuid=s.db_uuid AND i.session=s.session AND i.bundle=s.bundle AND i.level=s.level AND i.seq_idx=s.seq_idx+',offset,")")
        if(seqIdx<seqLen){
          selectStr=paste0(selectStr," || '->' || ")
        }
      }
    }
    selectStr=paste0(selectStr," END AS labels ")
  }else{
    selectStr=paste0(selectStr," '' AS labels ")
  }
  queryStr=paste(selectStr,fromStr,whereStr,orderStr,sep = ' ')
  
  # convert samples to milliseconds using SQL:
  seglist=DBI::dbGetQuery(emuDBhandle$connection, paste0("SELECT \
                       labels,
                       CASE type WHEN 'EVENT' THEN \
                        CAST (sample_start AS REAL)/ CAST( sample_rate AS REAL) * 1000.0 \
                       ELSE \
                        (CAST (sample_start AS REAL) + 0.5 ) / CAST( sample_rate AS REAL) * 1000.0 \
                       END AS start, \
                       CASE type WHEN 'EVENT' THEN \
                         0.0
                       ELSE \
                        (CAST (sample_end AS REAL) + 1.5 ) / CAST( sample_rate AS REAL) * 1000.0 \
                       END AS end, \
                       session || ':' || bundle AS utts, \
                       db_uuid,session,bundle, start_item_id AS startItemID, end_item_id AS endItemID,level,type,sample_start AS sampleStart,sample_end AS sampleEnd,sample_rate AS sampleRate\
                      FROM (", queryStr, ")"))
  
  # set emusegs type attribute, default 'segment'
  slType='segment'
  if(nrow(seglist)>0){
    # set to event only if all rows are of type EVENT
    dTypes=unique(seglist$level)
    if(length(dTypes)==1){
      if(dTypes[1]=='EVENT'){
        slType='event'
      }
    }
  }
  segmentList=make.emuRsegs(dbName = emuDBhandle$dbName, seglist = seglist,query = "FROM REQUERY", type = slType)
  return(segmentList)
}

##################################
##################################
##################################
equal.emusegs<-function(seglist1,seglist2,compareAttributes=TRUE,tolerance=0.0,uttsPrefix2=''){
  if(!inherits(seglist1,"emusegs")){
    stop("seglist1 is not of class emusegs")
  }
  if(!inherits(seglist2,"emusegs")){
    stop("seglist2 is not of class emusegs")
  }
  if(tolerance<0){
    stop("tolerance must be greater or equal 0")
  }
  
  sl1RowCnt=nrow(seglist1)
  sl2RowCnt=nrow(seglist2)
  if(sl1RowCnt!=sl2RowCnt){
    return(FALSE)
  }
  if(compareAttributes){
    attrEq=((attr(seglist1,'query')==attr(seglist2,'query'))&attr(seglist1,'database')==attr(seglist2,'database'))
    if(!attrEq){
      return(FALSE)
    }
  }
  # seglist have no implicit order
  sl1Order=order('[[.data.frame'(seglist1,'utts'),'[[.data.frame'(seglist1,'start'),'[[.data.frame'(seglist1,'end'),'[[.data.frame'(seglist1,'labels'))
  sl1=`[.data.frame`(seglist1,sl1Order,)
  sl2Order=order('[[.data.frame'(seglist2,'utts'),'[[.data.frame'(seglist2,'start'),'[[.data.frame'(seglist2,'end'),'[[.data.frame'(seglist2,'labels'))
  sl2=`[.data.frame`(seglist2,sl2Order,)
  equal=TRUE
  rs=1:sl1RowCnt
  for(i in rs){
    l1='[[.data.frame'(sl1,i,'labels')
    l2='[[.data.frame'(sl2,i,'labels')
    s1='[[.data.frame'(sl1,i,'start')
    s2='[[.data.frame'(sl2,i,'start')
    e1='[[.data.frame'(sl1,i,'end')
    e2='[[.data.frame'(sl2,i,'end')
    u1='[[.data.frame'(sl1,i,'utts')
    u2='[[.data.frame'(sl2,i,'utts')
    u2Sess=paste0(uttsPrefix2,u2)
    if(l1!=l2 | u1!=u2Sess){
      return(FALSE)
    }
    sdAbs=abs(s2-s1)
    if(sdAbs>tolerance){
      equal=FALSE
      #cat("Start differs ",s1,s2,sdAbs,"\n")
    }
    edAbs=abs(e2-e1)
    if(edAbs>tolerance){
      equal=FALSE
      #cat("End differs ",e1,e2,edAbs,"\n")
    }
    
  }
  return(equal)
  
}