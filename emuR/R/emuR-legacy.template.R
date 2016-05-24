requireNamespace("uuid", quietly = T)

## Create emuDB database schema object from EMU template (.tpl) file
## 
## @param tplPath EMU template file path
## @param dbUUID optional database UUID
## @param encoding  encoding of the template file
## @return object of class emuDB.schema.db
## @import stringr uuid wrassp
## @keywords emuDB database schema Emu 
## 
load_dbConfigFromEmuTemplate=function(tplPath,dbUUID=NULL,encoding=NULL){
  LEVEL_CMD='level'
  LABFILE_CMD='labfile'
  LABEL_CMD='label'
  SET_CMD='set'
  TRACK_CMD='track'
  PATH_CMD='path'
  LEGAL_CMD='legal'
  if(is.null(tplPath)) {
    stop("Argument tplPath (path to Emu template file) must not be NULL\n")
  }
  tplBasename = basename(tplPath)
  dbName=gsub("[.][tT][pP][lL]$","",tplBasename)
  
  # read
  if(is.null(encoding)){
    tpl = try(readLines(tplPath))
  }else{
    tpl = try(readLines(tplPath,encoding=encoding))
  }
  if(class(tpl) == "try-error") {
    stop("read tpl: cannot read from file ", tplPath)
  }
  # check if file (not directory)
  tplFInfo = try(file.info(tplPath))
  if(class(tplFInfo) == "try-error" | is.null(tplFInfo)) {
    stop("check template file: cannot get file info: ", tplPath)
  }
  if(tplFInfo[['isdir']]){
    stop(tplPath," is a directory. Expected a legacy EMU template file path.")
  }
  
  tracks=list()
  flags=list()
  levelDefinitions=list()
  linkDefinitions=list()
  pathDescriptors=list()
  annotationDescriptors=list()
  hlbTierDescriptors=list()
  hlbAnnotationDescriptor=NULL;
  
  lineNr=0
  
  for(line in tpl){
    lineNr=lineNr+1L
    trimmedLine=stringr::str_trim(line)
    if(trimmedLine!=''){
      firstChar=substr(trimmedLine,1,1)
      if(firstChar!='!'){
        
        lineTokensLst=strsplit(trimmedLine,'[[:space:]]+')
        lineTokens=lineTokensLst[[1]]
        lineTokenCount=length(lineTokens)
        if(lineTokenCount>=1){
          command=lineTokens[1]
          
          if(command==LABFILE_CMD){
            tierName=lineTokens[2]
            # TODO are there any default values for this properties?
            extension=NULL
            type=NULL
            timeFactor=NULL
            for(tki in 3:length(lineTokens)){
              tk=lineTokens[[tki]]
              
              if(substr(tk,1,1)==':'){
                # property key
                key=substring(tk,2)
              }else{
                #property value
                if(is.null(key)){
                  stop("Emu template parser key/value error in ",lineNr,"\n")
                }
                val=tk
                if(key=='extension'){
                  extension=val
                }else if(key=='type'){
                  type=val
                }else if(key=='time-factor'){
                  timeFactor=val
                }
                # reset key
                key=NULL 
              }
            }
            ad=list(name=tierName,extension=extension,type=type,timeFactor=timeFactor)
            annotationDescriptors[[length(annotationDescriptors)+1L]] <- ad
            
            # lab file can reference hlb level
            replaced=FALSE
            tdLen=length(levelDefinitions)
            for(i in 1:tdLen){
              td=levelDefinitions[[i]]
              if(td[['name']]==tierName){
                # replace 
                levelDefinitions[[i]]=list(name=td[['name']],type=type,attributeDefinitions=td[['attributeDefinitions']]);
                replaced=TRUE
                break;
              }
            }
            if(!replaced){
              # append
              levelDefinitions[[length(levelDefinitions)+1L]]=list(name=tierName);
            }
          }else if(command==TRACK_CMD){
            
            name=lineTokens[2]
            extension=lineTokens[3]
            track=list(name=name, columnName=name, fileExtension=extension)
            tracks[[length(tracks)+1L]] <- track
          }else if(command==SET_CMD){
            key=lineTokens[2]
            value=lineTokens[3]
            flags[[key]]=value
          }else if(command==PATH_CMD){
            annoKey=lineTokens[2]
            annoBasePath=lineTokens[3]
            pathDescr=list(basePath=annoBasePath,key=annoKey)
            pathDescriptors[[length(pathDescriptors)+1L]] <- pathDescr
            if(annoKey=='hlb'){
              # special meaning
              # hlb files are neither declared by tracks nor by labfile directive
              # add as annotationDescriptor
              ad=list(name=NULL,extension=annoKey,type='HLB')
              annotationDescriptors[[length(annotationDescriptors)+1L]] <- ad
            }
          }else if(command==LEVEL_CMD){
            levelTierName=lineTokens[2]
            
            if(lineTokenCount>=3){
              linkType="ONE_TO_MANY"
              if(lineTokenCount>=4){
                relationshipType=lineTokens[4]
                if(relationshipType=='many-to-many'){
                  linkType="MANY_TO_MANY"
                }
              }
              linkDefinition=list(type=linkType,superlevelName=lineTokens[3],sublevelName=levelTierName)
              linkDefinitions[[length(linkDefinitions)+1L]]=linkDefinition
            }
            tierDescr=list(name=levelTierName,type='ITEM', attributeDefinitions=list(list(name = levelTierName, type = "STRING")))
            exists=FALSE
            for(lDef in levelDefinitions){
              if(lDef[['name']]==levelTierName){
                exists=TRUE
                break
              }
            }
            if(!exists){
              levelDefinitions[[length(levelDefinitions)+1L]]=tierDescr
            }
            
            # TODO constraints
          }else if(command==LABEL_CMD){
            
            levelTierName=lineTokens[2]
            labelNames=list(levelTierName)
            if(lineTokenCount!=3){
              stop("Expected label directive \"label levelName labelName\"")
            }
            
            for(i in 1:length(levelDefinitions)){
              td=levelDefinitions[[i]]
              if(td[['name']]==levelTierName){
                # replace
                attrDefs=levelDefinitions[[i]][['attributeDefinitions']]
                attrDefs[[length(attrDefs)+1L]]=list(name=lineTokens[3],type='STRING')
                levelDefinitions[[i]]=list(name=levelTierName,type=td[['type']],attributeDefinitions=attrDefs);
                break
              }
            }
          }else if(command==LEGAL_CMD){
            if(lineTokenCount<=3){
              stop("Expected legal directive \"legal levelName groupName label1 label2 ... labeln\"")
            }
            attrName=lineTokens[2]
            labelGroupName=lineTokens[3]
            
            groupLabels=list()
            for(i in 4:lineTokenCount){
              groupLabels[[length(groupLabels)+1]]=lineTokens[i]
            }
            set=FALSE
            for(i in 1:length(levelDefinitions)){
              td=levelDefinitions[[i]]
              ads=td[['attributeDefinitions']]
              for(j in 1:length(ads)){
                ad=ads[[j]]
                if(ad[['name']]==attrName){
                  lblGrIdx=length(ad[['labelGroups']])+1
                  levelDefinitions[[i]][['attributeDefinitions']][[j]][['labelGroups']][[lblGrIdx]]=list(name=labelGroupName,values=groupLabels)
                  set=TRUE
                  break
                }
              }
              if(set){
                break
              }
            }
          }
        }
      }
    }
  }
  
  #pef=flags$PrimaryExtension
  tl=length(tracks)
  al=length(annotationDescriptors)
  # apply pathes to tracks  
  tss2=1:tl
  for(ti2 in tss2){
    for(pd in pathDescriptors){
      if(tracks[[ti2]][['fileExtension']] == pd[['key']]){
        tracks[[ti2]][['basePath']]=pd[['basePath']]
        break
      }
    }
  }
  
  # apply pathes to annotations
  as=1:al
  for(ai in as){
    for(pd in pathDescriptors){
      if(annotationDescriptors[[ai]][['extension']] == pd[['key']]){
        annotationDescriptors[[ai]][['basePath']]=pd[['basePath']]
        break
      }
    }
  }
  
  ssffTrackDefinitions=list()
  assign=list()
  mediafileBasePathPattern=NULL
  mediafileExtension=NULL
  for(tr in tracks){
    n=tr[['name']]
    e=tr[['fileExtension']]
    if(e==flags[['PrimaryExtension']]){
      primaryBasePath=tr[['basePath']]
    }
    if(n=='samples'){
      if(e!='wav'){
        cat("WARNING! Media file type with extension ",e," not supported by EMU-Webapp.\n")
      }
      mediafileExtension=e
      mediafileBasePathPattern=tr[['basePath']]
    }else{
      #array !
      ssffTrackDefinitions[[length(ssffTrackDefinitions)+1L]]=tr
      # default assign all to spectrum TODO
      
    }
  }
  
  if(is.null(dbUUID)){
  # Generate UUID 
  # problem: the UUID will change on every reload
  dbUUID=uuid::UUIDgenerate()
  }
  
  # default perspective
  # assign all SSFF tracks to sonagram
  assign=list()
  for(ssffTrack in ssffTrackDefinitions){
    # TODO dirty workaround
    # detect formant tracks by number of channels
    if(ssffTrack[['name']] == 'fm'){
      #ssffTrack$name='FORMANTS'
      #assign[[length(assign)+1]]=list(signalCanvasName='SPEC',ssffTrackName='FORMANTS')
    }
  }
  
  contourLims=list()
  sc=list(order=c("OSCI","SPEC"), assign=assign, contourLims=contourLims)
  
  defaultLvlOrder=list()
  for(ld in levelDefinitions){
    
    if(ld[['type']]=='SEGMENT' || ld[['type']]=='EVENT'){
      defaultLvlOrder[[length(defaultLvlOrder)+1L]]=ld[['name']]
    }
  }
  
  defPersp=list(name='default',signalCanvases=sc,levelCanvases=list(order=defaultLvlOrder),twoDimCanvases=list(order=list()))
  waCfg=list(perspectives=list(defPersp))
  dbSchema=list(name=dbName,UUID=dbUUID,mediafileBasePathPattern=mediafileBasePathPattern,mediafileExtension=mediafileExtension,ssffTrackDefinitions=ssffTrackDefinitions,levelDefinitions=levelDefinitions,linkDefinitions=linkDefinitions,EMUwebAppConfig=waCfg,annotationDescriptors=annotationDescriptors,tracks=tracks,flags=flags);
  
  # get max label array size
  maxLbls=0
  for(lvlDef in levelDefinitions){
    attrCnt=length(lvlDef[['attributeDefinitions']])
    if(attrCnt > maxLbls){
      maxLbls=attrCnt
    }
  }
  dbSchema[['maxNumberOfLabels']]=maxLbls
  return(dbSchema)
}