requireNamespace("stringr", quietly = T)

database.DDL.emuRsegsTmp = 'CREATE TEMP TABLE emursegs_tmp (
  labels TEXT,
  start FLOAT,
  end FLOAT,
  utts TEXT,
  db_uuid VARCHAR(36) NOT NULL,
  session TEXT, 
  bundle TEXT,
  start_item_id INTEGER,
  end_item_id INTEGER,
  level TEXT,
  type TEXT,
  sampleStart INTEGER,
  sampleEnd  INTEGER,
  sampleRate FLOAT
);'

create_requeryTmpTables <- function(emuDBhandle){
  DBI::dbGetQuery(emuDBhandle$connection, database.DDL.emuRsegsTmp)
}

drop_requeryTmpTables <- function(emuDBhandle){
  if("emursegs_tmp" %in% DBI::dbListTables(emuDBhandle$connection)){
    DBI::dbGetQuery(emuDBhandle$connection, "DROP TABLE emursegs_tmp")
  }
}


##' Requery sequential context of segment list in an emuDB
##' @description Function to requery sequential context of a segment list queried from an emuDB
##' @details Builds a new segment list on the same hierarchical level and the same length as the segment list given in \code{seglist}. The resulting segments usually have different start position and length (in terms of items of the respective level) controlled by the \code{offset},\code{offsetRef} and \code{length} parameters.
##' A segment here is defined as a single item or a chain of items from the respective level, e.g. if a level in a bundle instance has labels 'a', 'b' and 'c' in that order, 'a' or 'a->b' oder 'a->b->c' are all valid segments, but not 'a->c'.
##' \code{offsetRef} determines if the position offset is referenced to the start or the end item of the segments in the input list \code{seglist}; parameter \code{offset} determines the offset of the resulting item start position to this reference item;
##' parameter \code{length} sets the item length of the result segments.
##' If the requested segments are out of bundle item boundaries and parameter \code{ignoreOutOfBounds} is \code{FALSE} (the default), an error is generated. To get residual resulting segments that lie within the bounds the \code{ignoreOutOfBounds} parameter can be set to \code{TRUE}.
##' The returned segment list is usually of the same length and order as the input \code{seglist}; if \code{ignoreOutOfBounds=FALSE}, the resulting segment list may be out of sync.
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param seglist segment list to requery on (type: 'emuRsegs')
##' @param offset start item offset in sequence (default is 0, meaning the start or end item of the input segment)
##' @param offsetRef reference item for offset: 'START' for first and 'END' for last item of segment
##' @param length item length of segments in the returned segment list
##' @param ignoreOutOfBounds ignore result segments that are out of bundle bounds
##' @return result set object of class 'emuRsegs' containing the requeried segments
##' @export
##' @seealso \code{\link{query}} \code{\link{requery_hier}} \code{\link{emuRsegs}}
##' @keywords emuDB database requery
##' @examples
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' ## Requery previous item of 'p' on level 'Phonetic'
##' sl1 = query(ae, "Phonetic == p")
##' 
##' requery_seq(ae, sl1, offset = -1)
##' 
##' ## Requery context (adding previuos and following elements) of 'p' on phonetic level
##'
##' requery_seq(ae, sl1, offset = -1, length = 3)
##' 
##' ## Requery previous item of n->t sequence
##' sl2 = query(ae, "[Phoneme == n -> Phoneme == t]")
##' 
##' requery_seq(ae, sl2, offset = -1)
##' 
##' ## Requery last item within n->t sequence
##' 
##' requery_seq(ae, sl2, offsetRef = 'END')
##' 
##' ## Requery following item after n->t sequence
##' 
##' requery_seq(ae, sl2, offset = 1, offsetRef = 'END')
##' 
##' ## Requery context (previous and following items) of n->t sequence
##' 
##' requery_seq(ae, sl2, offset = -1, length = 4)
##' 
##' ## Requery next word contexts (sequence includes target word)
##' 
##' sl3 = query(ae, "Text == to")
##' requery_seq(ae, sl3, length = 2)
##' 
##' ## Requery following two word contexts, ignoring segment 
##' ## sequences that are out of bundle end bounds 
##' requery_seq(ae, sl3, length = 3, ignoreOutOfBounds = TRUE)
##' 
##' }
requery_seq<-function(emuDBhandle, seglist, offset=0,offsetRef='START',length=1,ignoreOutOfBounds=FALSE){
  if(!inherits(seglist,"emuRsegs")){
    stop("Segment list 'seglist' must be of type 'emuRsegs'. (Do not set a value for 'resultType' parameter in the query() command; then the default resultType=emuRsegs will be used)")
  }
  if(length<=0){
    stop("Parameter length must be greater than 0")
  }

  if(nrow(seglist)==0){
    # empty seglist, return the empty list
    return(seglist)
  }else{
    # drop create tmp tables and recreate (will ensure they are empty)
    drop_requeryTmpTables(emuDBhandle)
    create_requeryTmpTables(emuDBhandle)
    # place in emuRsegsTmp table
    DBI::dbGetQuery(emuDBhandle$connection, "DELETE FROM emursegs_tmp;") # delete 
    DBI::dbWriteTable(emuDBhandle$connection, "emursegs_tmp", as.data.frame(seglist), append=T, row.names = F) # append to make sure field names done't get overwritten
    
    # load config
    dbConfig=load_DBconfig(emuDBhandle)
    
    # query for sequential requeries
    heQueryStr=paste0("SELECT il.db_uuid,il.session,il.bundle,il.item_id AS seq_start_id,ir.item_id AS seq_end_id,",length," AS seq_len,sl.level FROM emursegs_tmp sl,items sll,items slr,items il, items ir \
                        WHERE \
                         il.db_uuid=ir.db_uuid AND il.session=ir.session AND il.bundle=ir.bundle AND \
                         il.db_uuid=sl.db_uuid AND il.session=sl.session AND il.bundle=sl.bundle AND \
                         sll.db_uuid=sl.db_uuid AND sll.session=sl.session AND sll.bundle=sl.bundle AND sl.start_item_id=sll.item_id AND \
                         slr.db_uuid=sl.db_uuid AND slr.session=sl.session AND slr.bundle=sl.bundle AND sl.end_item_id=slr.item_id AND ")
    if(offsetRef=='START'){
      heQueryStr=paste0(heQueryStr,"il.level=sll.level AND il.seq_idx=sll.seq_idx+",offset," AND \
                          ir.level=sll.level AND ir.seq_idx=sll.seq_idx+",offset+length-1)
    }else if(offsetRef=='END'){
      heQueryStr=paste0(heQueryStr,"il.level=slr.level AND il.seq_idx=slr.seq_idx+",offset," AND \
                          ir.level=slr.level AND ir.seq_idx=slr.seq_idx+",offset+length-1)
    }else{
      stop("Parameter offsetRef must be one of 'START' or 'END'\n")
    }
    heQueryStr=paste0(heQueryStr," ORDER BY il.ROWID");
    he = DBI::dbGetQuery(emuDBhandle$connection, heQueryStr)
    slLen=nrow(seglist)
    resLen=nrow(he)
    outOfBndCnt=slLen-resLen
    if(!ignoreOutOfBounds & outOfBndCnt>0){
      if(outOfBndCnt==slLen){
        stop("All (",outOfBndCnt,") of the requested sequence(s) is/are out of boundaries.")
      }else{
        stop(outOfBndCnt," of the requested sequence(s) is/are out of boundaries.\nSet parameter 'ignoreOutOfBounds=TRUE' to get residual result segments that lie within the bounds.")
      }
    }
    
    # drop and create tmpQueryTables and write to table
    drop_allTmpTablesDBI(emuDBhandle)
    create_tmpFilteredQueryTablesDBI(emuDBhandle)
    DBI::dbWriteTable(emuDBhandle$connection, "interm_res_items_tmp_root", he, overwrite=T)
    

    trSl=convert_queryResultToVariableEmuRsegs(emuDBhandle)
    drop_requeryTmpTables(emuDBhandle)
    
    return(trSl)
  }
}

##' Requery hierarchical context of a segment list in an emuDB
##' @description Function to requery hierarchical context of a segment list queried from an emuDB
##' @details A segment is defined as a single item or a chain of items from the respective level, e.g. if a level in a bundle instance has labels 'a', 'b' and 'c' in that order, 'a' or 'a->b' oder 'a->b->c' are all valid segments, but not 'a->c'.
##' For each segment of the input segment list \code{seglist} the function checks the start and end item for hierarchically linked items in the given target level, and based on them constructs segments in the target level.
##' As the start item in the resulting segment the item with the lowest sample position is chosen; for the end item that with the highest sample position.
##' If result and input segment list have the same length (for each input segment one segment on the target level was found), the result segment list has the same length and order as the input list; 
##' in 'upwards' requeries this can cause a resulting segment list to contain two (or more) copies of the same segment, if the same item from the input list was linked twice or more to an item of the target level, e.g. a phoneme 'p' requeried to the word level might result in two identical segments 'Papa' in the result list. 
##' If the length of input and output list differ (e.g. because a link is missing in the emuDB), a synchronous ordering is not possible and therefore a warning is generated.
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param seglist segment list to requery on (type: \link{emuRsegs})
##' @param level character string: name of target level 
##' @return result set object of class \link{emuRsegs}
##' @export
##' @seealso \code{\link{query}} \code{\link{requery_seq}} \code{\link{emuRsegs}}
##' @keywords emuDB database requery
##' @examples
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' ## Downward requery: find 'Phoneme' sequences of all words 'beautiful' (of level 'Text')
##' ## Note that the resulting segments consists of phoneme sequences and have therefore 
##' ## the same length as the word segments.
##'
##' sl1 = query(ae, "Text == beautiful")
##' requery_hier(ae, sl1, level = "Phoneme")
##'
##' ## Upward requery: find all word segments that dominate a 'p' on level 'Phoneme'
##' ## Note that the resulting segments are larger than the input segments,
##' ## because they contain the complete words.
##' 
##' sl1 = query(ae, "Phonetic == p")
##' requery_hier(ae, sl1, level = 'Text')
##' 
##' ## Why is there a 'p' the word 'emphazised'? Requery the whole words back down to 'Phoneme' level:
##'
##' requery_hier(ae, sl1, level = 'Phoneme')
##'
##' ## ... because of 'stop epenthesis' a 'p' is inserted between 'm' and 'f'
##' 
##' ## Combined requery: last phonemes of all words beginning with 'an'.
##' ## Note that we use a regular expression 'an.*' (EQL operator '=~') in the query.
##' 
##' sl1=query(ae, "Text =~ an.*")
##' requery_seq(ae, requery_hier(ae, sl1, level = 'Phoneme'), offsetRef = 'END')
##' 
##' }
requery_hier<-function(emuDBhandle, seglist, level=NULL){
  if(!inherits(seglist,"emuRsegs")){
    stop("Segment list 'seglist' must be of type 'emuRsegs'. (Do not set a value for 'resultType' parameter for the query, the default resultType will be used)")
  }

  if(nrow(seglist)==0){
    # empty seglist, return the empty list
    return(seglist)
  }else{
    # drop create tmp tables and recreate (will ensure they are empty)
    drop_requeryTmpTables(emuDBhandle)
    create_requeryTmpTables(emuDBhandle)
    # place in emursegs_tmp table
    DBI::dbGetQuery(emuDBhandle$connection, "DELETE FROM emursegs_tmp;")
    DBI::dbWriteTable(emuDBhandle$connection, "emursegs_tmp", as.data.frame(seglist), append=T, row.names = F) # append to avoid rewirte of col names
    
    # load config
    dbConfig=load_DBconfig(emuDBhandle)
    
    targetRootLevelName=NULL
    if(is.null(level)){
      heQueryStr=paste0("SELECT il.db_uuid,il.session,il.bundle,il.item_id AS seq_start_id,ir.item_id AS seq_end_id,ir.seq_idx-il.seq_idx+1 AS seq_len,il.level \
                          FROM \
                          ( SELECT ils.*,min(ils.seq_idx),sl.ROWID AS lrId FROM items ils,items slil,seglist sl WHERE \
                          ils.db_uuid=sl.db_uuid AND ils.session=sl.session AND ils.bundle=sl.bundle AND \
                          slil.db_uuid=sl.db_uuid AND slil.session=sl.session AND slil.bundle=sl.bundle AND \
                          slil.item_id=sl.start_item_id AND ils.level=slil.level AND (\
                          (ils.item_id=sl.start_item_id) OR 
                          (EXISTS (SELECT * FROM links_ext lr \
                          WHERE lr.db_uuid=sl.db_uuid AND lr.session=sl.session AND lr.bundle=sl.bundle \
                          AND ((lr.from_id=sl.start_item_id AND lr.to_id=ils.item_id) OR (lr.from_id=ils.item_id AND lr.to_id= sl.start_item_id))\
                          )) \
                          ) GROUP BY lrId ) \
                          AS il JOIN \
                          ( SELECT irs.*,max(irs.seq_idx),sl.ROWID AS rrId FROM items irs,items slir,seglist sl WHERE \
                          irs.db_uuid=sl.db_uuid AND irs.session=sl.session AND irs.bundle=sl.bundle AND \
                          slir.db_uuid=sl.db_uuid AND slir.session=sl.session AND slir.bundle=sl.bundle AND \
                          slir.item_id=sl.end_item_id AND irs.level=slir.level AND (\
                          (irs.item_id=sl.end_item_id) OR
                          (EXISTS (SELECT * FROM links_ext lr \
                          WHERE lr.db_uuid=sl.db_uuid AND lr.session=sl.session AND lr.bundle=sl.bundle \
                          AND ((lr.from_id=sl.end_item_id AND lr.to_id=irs.item_id) OR (lr.from_id=irs.item_id AND lr.to_id= sl.end_item_id))\
                          )) \
                          ) GROUP BY rrId ) \
                          AS ir ON lrId=rrId
                          ")
      
    }else{
      
      check_levelAttributeName(emuDBhandle,level)
      targetRootLevelName=get_levelNameForAttributeName(emuDBhandle, attributeName = level)
      heQueryStr=paste0("SELECT il.db_uuid,il.session,il.bundle,il.item_id AS seq_start_id,ir.item_id AS seq_end_id,(ir.seq_idx-il.seq_idx+1) AS seq_len,'",level,"' AS level \
                              FROM 
                              ( SELECT ils.*,min(ils.seq_idx),sll.ROWID AS lrId FROM emursegs_tmp sll,items ils WHERE \
                              ils.db_uuid=sll.db_uuid AND ils.session=sll.session AND ils.bundle=sll.bundle AND \
                              ils.level='",targetRootLevelName,"' AND (\
                              (ils.item_id=sll.start_item_id) OR 
                              (EXISTS (SELECT * FROM links_ext ll \
                              WHERE ll.db_uuid=sll.db_uuid AND ll.session=sll.session AND ll.bundle=sll.bundle \
                                  AND ((ll.from_id=sll.start_item_id AND ll.to_id=ils.item_id) OR (ll.from_id=ils.item_id AND ll.to_id= sll.start_item_id))\
                                  )) \
                              ) GROUP BY lrId ORDER BY lrId,ils.seq_idx) \
                              AS il JOIN \
                              ( SELECT irs.*,max(irs.seq_idx),slr.ROWID AS rrId FROM emursegs_tmp slr,items irs WHERE \
                              irs.db_uuid=slr.db_uuid AND irs.session=slr.session AND irs.bundle=slr.bundle AND \
                              irs.level='",targetRootLevelName,"' AND (\
                              (irs.item_id=slr.end_item_id) OR
                              (EXISTS (SELECT * FROM links_ext lr \
                              WHERE lr.db_uuid=slr.db_uuid AND lr.session=slr.session AND lr.bundle=slr.bundle \
                                  AND ((lr.from_id=slr.end_item_id AND lr.to_id=irs.item_id) OR (lr.from_id=irs.item_id AND lr.to_id= slr.end_item_id))\
                                )) \
                              ) GROUP BY rrId ORDER BY rrId,irs.seq_idx DESC) \
                              AS ir ON lrId=rrId ")
      
    }

    he = DBI::dbGetQuery(emuDBhandle$connection, heQueryStr)
    
    # drop and create tmpQueryTables and write to table
    drop_allTmpTablesDBI(emuDBhandle)
    create_tmpFilteredQueryTablesDBI(emuDBhandle)
    DBI::dbWriteTable(emuDBhandle$connection, "interm_res_items_tmp_root", he, overwrite=T)
    
    trSl=convert_queryResultToVariableEmuRsegs(emuDBhandle)
    inSlLen=nrow(seglist)
    trSlLen=nrow(trSl)

    if(inSlLen!=trSlLen){
      warning("Length of requery segment list (",trSlLen,") differs from input list (",inSlLen,")!")
    }
    return(trSl)
  }
}
