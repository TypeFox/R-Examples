requireNamespace("stringr", quietly = T)
requireNamespace("RSQLite", quietly = T)

## Parser for Praat TextGrid files
## 
## parses directly to DBI tables (items, labels)
## @param emuDBhandle 
## @param textGridPath TextGrid file connection
## @param sampleRate sample rate of correponding signal file
## @param encoding text encoding (currently the only excepted is the default UTF-8)
## @param bundle name of bundle 
## @param session name of session
## 
parse_TextGridDBI <- function(emuDBhandle, TextGridPath=NULL, sampleRate, encoding="UTF-8", 
                              bundle=NULL, session="0000") {
  
  #####################
  # check arguments (TODO better checks for classes and the like...)
  
  if(is.null(TextGridPath)) {
    stop("Argument TextGridPath must not be NULL\n")
  }
  if(sampleRate <=0 ){
    stop("Samplerate must be greater than zero\n")
  }
  if(encoding != "UTF-8"){
    stop("The only encoding that is currently supported is UTF-8\n")
  }
  if(is.null(bundle)){
    stop("Argument bundle must not be NULL!\n")
  }
  if(is.null(session)){
    stop("Argument session must not be NULL!\n")
  }
  
  #
  #####################
  
  itemCounterGlobal = 1
  itemCounterLevel = 1
  
  FILE_TYPE_KEY="File type"
  OBJECT_CLASS_KEY="Object class"
  TIERS_SIZE_KEY="size"
  TIER_ITEM_KEY="item"
  NAME_KEY="name"
  INTERVALS_KEY="intervals"
  POINTS_KEY="points"
  XMIN_KEY="xmin"
  XMAX_KEY="xmax"
  TEXT_KEY="text"
  TIME_KEY="time"
  
  FILE_TYPE_VAL_OO_TEXTFILE="ooTextFile"
  OBJECT_CLASS_VAL_TEXTGRID="TextGrid"
  TIER_CLASS_VAL_INTERVAL="IntervalTier"
  TIER_CLASS_VAL_TEXT="TextTier"
  
  fileType=NULL
  objectClass=NULL
  hasTiers=FALSE
  tiersCount=NULL
  currentTier=NULL
  currentTierClass=NULL
  currentTierName=NULL
  currentTierSize=NULL
  
  # read TextGrid
  tg = try(readLines(TextGridPath))
  if(class(tg) == "try-error") {
    stop("read.TextGrid: cannot read from file ", TextGridPath)
  }
  
  # remove all trailing/leading white spaces (for speed improvment)
  tg = gsub("^\\s+|\\s+$", "", tg)
  
  for(line in tg){
    # check for fileType
    if(is.null(fileType)){
      p=parse_lineToKeyValue(line,doubleQuoted=TRUE, initialTrim=FALSE)
      if(! is.null(p)){
        if(p[1]==FILE_TYPE_KEY){
          fileType=p[2]
          # check if of correct type:
          if(fileType != FILE_TYPE_VAL_OO_TEXTFILE){
            stop("Can only parse TextGrids with the File type: ", FILE_TYPE_VAL_OO_TEXTFILE, 
                 ". Found following File type: ", fileType)
          }
          
        }
      }
    }else{
      # check for objectClass
      if(is.null(objectClass)){
        p=parse_lineToKeyValue(line,doubleQuoted=TRUE, initialTrim=FALSE)
        if(! is.null(p)){
          if(p[1]==OBJECT_CLASS_KEY){
            objectClass=p[2]
          }
        }
      }else{
        # if we have both the file type and the object class        
        if((fileType==FILE_TYPE_VAL_OO_TEXTFILE) && (objectClass==OBJECT_CLASS_VAL_TEXTGRID)){
          
          if(is.null(tiersCount)){
            
            p=parse_lineToKeyValue(line, initialTrim=FALSE)
            if((!is.null(p)) && (p[1]=='size')){
              tiersCount=p[2]
            }
          }else{
            ## if we have tiersCount tiers
            if(length(grep("^item",line))==1){
              
              tierIndexStr=sub('item\\s*','',line);
              tierIndexStr=sub('\\s*:$','',tierIndexStr);
              if(length(grep('\\[\\s*[0-9]+\\s*\\]',tierIndexStr))==1){
                tierIndexStr=sub('\\[\\s*','',tierIndexStr);
                tierIndexStr=sub('\\s*\\]','',tierIndexStr);
                
                tierIndex=tierIndexStr;
                # reset level/tier attributes
                itemCounterLevel = 1
                currentTierClass=NULL;
                currentTierName=NULL;
                currentTierSize=NULL;
                currentSegment=NULL;
                currentSegmentIndex=NULL;
                currentSegmentStart=NULL;
                currentSegmentDur=NULL;
                currentSegmentLabel=NULL;
                currentMark=NULL;
                currentPointIndex=NULL;
                currentPointSample=NULL;
                currentPointLabel=NULL;
              }
            }else {
              # check for currentTierClass
              if(is.null(currentTierClass)){
                p=parse_lineToKeyValue(line,doubleQuoted=TRUE, initialTrim=FALSE)
                if((! is.null(p)) && ('class' == p[1])){
                  currentTierClass=p[2];
                  if(currentTierClass==TIER_CLASS_VAL_INTERVAL){
                    
                  }else if(currentTierClass==TIER_CLASS_VAL_TEXT){
                  }else{
                    stop("TextGrid tiers of class \"",currentTierClass,"\" not supported!");
                  }
                }
              }
              # check for currentTierName
              if(is.null(currentTierName)){
                p=parse_lineToKeyValue(line,doubleQuoted=TRUE, initialTrim=FALSE)
                if((! is.null(p)) && ('name' == p[1])){
                  currentTierName=p[2];
                  
                }
              }
              # if we have the currentTierClass
              if(!is.null(currentTierClass)){
                if(currentTierClass==TIER_CLASS_VAL_INTERVAL){
                  # find size (and other properties)
                  if((is.null(currentTierSize)) && (length(grep('^intervals[[:space:]]*:.*',line))==1)){
                    
                    intervalsPropertyStr=stringr::str_trim(sub('^intervals[[:space:]]*:','',line))
                    intervalsProperty=parse_lineToKeyValue(intervalsPropertyStr,initialTrim=FALSE);
                    if((!is.null(intervalsProperty)) && (intervalsProperty[1]=='size')){
                      currentTierSize=intervalsProperty[2]
                      #cat("intervals: size=",currentTierSize,"\n");
                      
                    }
                  }
                  if(length(grep('intervals[[:space:]]*[[][[:space:]]*[0-9]+[[:space:]]*[]][[:space:]]*[:][[:space:]]*',line))==1){
                    
                    segmentIndexStr=sub("intervals[[:space:]]*[[][[:space:]]*","",line);
                    segmentIndexStr=sub("[[:space:]]*[]][[:space:]]*[:][[:space:]]*","",segmentIndexStr);
                    currentElementIndex=segmentIndexStr;
                    
                    currentSegmentIndex=segmentIndexStr;
                    currentSegmentStart=NULL;
                    currentSegmentEnd=NULL;
                    currentSegmentLabel=NULL;
                  }else{
                    p=parse_lineToKeyValue(line, doubleQuoted=TRUE, initialTrim=FALSE)
                    if((!is.null(p)) && (!is.null(currentSegmentIndex))){
                      if(p[1] == "xmin"){
                        minTimeStr=p[2]
                        minTime=as(minTimeStr,"numeric")
                        startSample = floor(minTime * sampleRate)
                        currentSegmentStart=startSample
                      }else if(p[1]=="xmax"){
                        maxTimeStr=p[2];
                        maxTime=as(maxTimeStr,"numeric")
                        currentSegmentEnd = floor(maxTime * sampleRate)
                        
                        
                      }else if(p[1]=="text"){
                        label=p[2];
                        currentSegmentLabel=label;
                      }
                      
                      if(!is.null(currentSegmentIndex) && 
                         !is.null(currentSegmentStart) &&
                         !is.null(currentSegmentEnd) &&
                         !is.null(currentSegmentLabel)){
                        sampleDur = currentSegmentEnd - currentSegmentStart - 1
                        labels=list(list(name=currentTierName,value=currentSegmentLabel))
                        
                        
                        # item entry:
                        DBI::dbGetQuery(emuDBhandle$connection, paste0("INSERT INTO items VALUES"," ('", emuDBhandle$UUID, "', '", session, "', '", bundle, "', '", itemCounterGlobal, 
                                                                   "', '", currentTierName, "', '", "SEGMENT", 
                                                                   "', ", itemCounterLevel, ", ", sampleRate, ", ", "NULL", ", ", currentSegmentStart, 
                                                                   ", ", sampleDur, ")"))
                        
                        
                        
                        # label entry:
                        DBI::dbGetQuery(emuDBhandle$connection, paste0("INSERT INTO labels VALUES","('", 
                                                                   emuDBhandle$UUID, "', '", session, "', '", bundle, "',", itemCounterGlobal,
                                                                   ", ", 0,", '", currentTierName, "', '", gsub("'","''", currentSegmentLabel), "')"))
                        
                        # links entry:
                        # no link entry because TextGrids don't have hierarchical infos
                        
                        # increase counters
                        itemCounterGlobal = itemCounterGlobal + 1
                        itemCounterLevel = itemCounterLevel + 1
                        
                        
                        currentSegment=NULL;
                        currentSegmentIndex=NULL;
                        currentSegmentStart=NULL;
                        currentSegmentDur=NULL;
                      }
                      
                    }
                  }
                  
                }else if(currentTierClass==TIER_CLASS_VAL_TEXT){
                  # find size (and other properties)
                  if((is.null(currentTierSize)) && (length(grep('^points[[:space:]]*[:].*',line))==1)){
                    
                    intervalsPropertyStr=stringr::str_trim(sub('^points[[:space:]]*[:]','',line))
                    intervalsProperty=parse_lineToKeyValue(intervalsPropertyStr, initialTrim=FALSE);
                    if((!is.null(intervalsProperty)) && (intervalsProperty[1]=='size')){
                      currentTierSize=intervalsProperty[2]
                    }
                  }
                  if(length(grep("points[[:space:]]*[[][[:space:]]*[0-9]+[[:space:]]*[]][[:space:]]*[:][[:space:]]*",line))==1){
                    pointIndexStr=sub("points[[:space:]]*[[][[:space:]]*","",line);
                    pointIndexStr=sub("[[:space:]]*[]][[:space:]]*[:][[:space:]]*","",pointIndexStr);
                    currentPointIndex=as.integer(pointIndexStr)
                    currentElementIndex=currentPointIndex
                    currentPointLabel=NULL;
                    currentPointSample=NULL;
                  }else{
                    #cat("inside point: \n")
                    p=parse_lineToKeyValue(line,doubleQuoted=TRUE, initialTrim=FALSE)
                    if((!is.null(p)) && (!is.null(currentPointIndex))){
                      if(p[1]=="time" || p[1]=="number"){
                        timePointStr=p[2];
                        timePoint=as(timePointStr,"numeric")
                        samplePoint = floor(timePoint * sampleRate)
                        currentPointSample=samplePoint
                      }else if(p[1]=="mark"){
                        currentPointLabel=p[2]
                      }else if(p[1]=="text"){
                        currentPointLabel=p[2]
                      }
                    }
                    if(!is.null(currentPointIndex) && 
                       !is.null(currentPointSample) &&
                       !is.null(currentPointLabel)){
                      
                      labels=list(list(name=currentTierName,value=currentPointLabel))
                      
                      # item entry
                      itemId = paste0(emuDBhandle$dbName, '_', session, '_', bundle, '_', itemCounterGlobal)
                      
                      
                      DBI::dbGetQuery(emuDBhandle$connection, paste0("INSERT INTO items VALUES"," ('", emuDBhandle$UUID, "', '", session, "', '", bundle, "', ",
                                                                 itemCounterGlobal, ", '", currentTierName,"', '", "EVENT", 
                                                                 "', ", itemCounterLevel, ", ", sampleRate, ", ", currentPointSample, ", ", "NULL", 
                                                                 ", ", "NULL", ")"))
                      
                      
                      # label entry:
                      DBI::dbGetQuery(emuDBhandle$connection, paste0("INSERT INTO labels VALUES","('", 
                                                                 emuDBhandle$UUID, "', '", session, "', '", bundle, "',", itemCounterGlobal,
                                                                 ", ", 0,", '", currentTierName, "', '", gsub("'","''", currentPointLabel), "')"))              
                      
                      
                      
                      # links entry:
                      # no link entry because TextGrids don't have hierarchical infos
                      
                      # increase counters
                      itemCounterGlobal = itemCounterGlobal + 1
                      itemCounterLevel = itemCounterLevel + 1
                      
                      currentPointIndex=NULL;
                      currentPointLabel=NULL;
                      currentPointSample=NULL;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

#############################
TextGridToBundleAnnotDFs <- function(tgPath, sampleRate, name, annotates){
  
  FILE_TYPE_KEY="File type"
  OBJECT_CLASS_KEY="Object class"
  
  tgChar = readChar(tgPath, file.info(tgPath)$size)
  lines = unlist(strsplit(tgChar, "\n"))
  
  if(!grepl(paste0("^", FILE_TYPE_KEY), lines[1]) & !grepl(paste0("^", OBJECT_CLASS_KEY), lines[2])){
    stop("First two lines of TextGrid file do not match: ", FILE_TYPE_KEY, "; and: ", OBJECT_CLASS_KEY, ". Only long form TextGrids are currently supported. Problem file is: ", tgPath)
  }
  
  
  # estimate how many items are in TextGrid for preallocation
  nrOfItems = length(grep("text|mark\\s*=", lines))
  
  # init data frames and preallocate enough rows
  items = data.frame(item_id = integer(nrOfItems), 
                     level = character(nrOfItems),
                     type = character(nrOfItems),
                     seq_idx = integer(nrOfItems),
                     sample_rate = numeric(nrOfItems),
                     sample_point = integer(nrOfItems),
                     sample_start = integer(nrOfItems),
                     sample_dur = integer(nrOfItems),
                     stringsAsFactors = F)
  
  labels = data.frame(item_id = integer(nrOfItems),
                      label_idx = integer(nrOfItems),
                      name = character(nrOfItems),
                      label = character(nrOfItems),
                      stringsAsFactors = F)
  
  
  
  # split at "...items [1]..." type lines
  tiers = unlist(strsplit(tgChar, ".*item\\s\\[[0-9]+\\].*\n", perl = T))
  header = tiers[1] # extract header
  tiers = tiers[-1]
  
  maxItemID = 1
  # iterate through tiers
  for(i in 1:length(tiers)){
    curTier = tiers[i]
    tierLines = unlist(strsplit(curTier, "\n"))
    
    tierHeaderEndIdx = grep("[intervals|points]:\\s*size", tierLines, perl = T)
    tierHeader = tierLines[1:tierHeaderEndIdx]
    tierLines = tierLines[-1:(-1*tierHeaderEndIdx)]
    
    
    if(grepl("IntervalTier", tierHeader[1])){
      
      levelName = sub('\\"\\s*$', "", sub('^\\s*name\\s*=\\s*\\"', "", tierHeader[grepl("name", tierHeader)], perl = T), perl = T)
      xminTimes = as.numeric(sub("^\\s*xmin\\s*=\\s*", "", tierLines[grepl("xmin", tierLines)], perl = T)) # as.numeric seems to be able to deal with trailing blanks
      xmaxTimes = as.numeric(sub("^\\s*xmax\\s*=\\s*", "", tierLines[grepl("xmax", tierLines)], perl = T)) # as.numeric seems to be able to deal with trailing blanks
      texts = sub('\\"\\s*$', "", sub('^\\s*text\\s*=\\s*\\"', "", tierLines[grepl("text", tierLines)]), perl = T)
      
      # calculate times 
      startSamples = floor(xminTimes * sampleRate)
      endSamples = floor(xmaxTimes * sampleRate)
      sampleDurs = endSamples - startSamples - 1
      
      # insert in data frames  
      items[maxItemID:(maxItemID + length(xminTimes) - 1), ]  = data.frame(item_id = maxItemID:(maxItemID + length(xminTimes) - 1), 
                                                                           level = rep(levelName, length(xminTimes)),
                                                                           type = rep("SEGMENT", length(xminTimes)),
                                                                           seq_idx = 1:length(xminTimes),
                                                                           sample_rate = rep(sampleRate, length(xminTimes)),
                                                                           sample_point = NA,
                                                                           sample_start = startSamples,
                                                                           sample_dur = sampleDurs,
                                                                           stringsAsFactors = F)
      
      labels[maxItemID:(maxItemID + length(xminTimes) - 1), ] = data.frame(item_id = maxItemID:(maxItemID + length(xminTimes) - 1),
                                                                           label_idx = rep(1, length(xminTimes)),
                                                                           name = rep(levelName, length(xminTimes)),
                                                                           label = texts,
                                                                           stringsAsFactors = F)
      
      maxItemID = max(items$item_id) + 1
      
    }else if(grepl("TextTier", tierHeader[1])){
      levelName = sub('\\"\\s*$', "", sub('^\\s*name\\s*=\\s*\\"', "", tierHeader[grepl("name", tierHeader)], perl = T), perl = T)
      pointsTimes = as.numeric(sub("^\\s*number\\s*=\\s*", "", tierLines[grepl("number", tierLines)], perl = T)) # as.numeric seems to be able to deal with trailing blanks
      marks = sub('\\"\\s*$', "", sub('^\\s*mark\\s*=\\s*\\"', "", tierLines[grepl("mark", tierLines)]), perl = T)
      
      # calculate times 
      samplePoints = floor(pointsTimes * sampleRate)
      
      # create data frames
      items[maxItemID:(maxItemID + length(samplePoints) - 1), ] = data.frame(item_id = maxItemID:(maxItemID + length(pointsTimes) - 1), 
                         level = rep(levelName, length(pointsTimes)),
                         type = rep("EVENT", length(pointsTimes)),
                         seq_idx = 1:length(pointsTimes),
                         sample_rate = rep(sampleRate, length(pointsTimes)),
                         sample_point = samplePoints,
                         sample_start = NA,
                         sample_dur = NA,
                         stringsAsFactors = F)
      
      labels[maxItemID:(maxItemID + length(samplePoints) - 1), ] = data.frame(items_id = maxItemID:(maxItemID + length(pointsTimes) - 1),
                          label_idx = rep(1, length(pointsTimes)),
                          name = rep(levelName, length(pointsTimes)),
                          label = marks,
                          stringsAsFactors = F)
      
      maxItemID = max(items$item_id) + 1
      
    }else{
      stop("Found Tier that does not have a class definition 'IntervalTier' or 'TextTier'. This probably means it is a mal formated TextGrid file. Problem file is: ", tgPath)
    }
    
  }
  
  links = data.frame(bundle = character(), from_id = integer(), to_id = integer(), label = character(), stringsAsFactors = F)
  
  return(list(name = name, annotates = annotates, sampleRate = sampleRate, items = items, links = links, labels = labels))
}

# FOR DEVELOPMENT
# library('testthat')
# test_file('tests/testthat/test_aaa_initData.R')
# test_file('tests/testthat/test_emuR-parse_TextGrid.R')
# tgPath = "~/Desktop/emuR_demoData/TextGrid_collection/msajc003.TextGrid"

