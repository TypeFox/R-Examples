build_legacyBundleList<-function(parsedEmuPath,currentPath=NULL,fileSuffixPattern,bundleList=list()){
  if(length(parsedEmuPath)==0){
    fileRegexPattern=paste0('^.*',fileSuffixPattern)
    fileList = list.files(currentPath, pattern=fileRegexPattern, recursive=F, full.names=F)
    if(length(fileList)==0){
      return(NULL)
    }else{
      # TODO convert to bundle names
      res=list()
      for(f in fileList){
        # remove file suffix for bundle name
        bundleName=gsub(x = f,pattern=fileSuffixPattern,replacement='')
        res[[length(res)+1]]=bundleName
      }
      return(res)
    }
  }else{
    p1=parsedEmuPath[[1]]
    if(length(parsedEmuPath)==1){
      restPath=list()
    }else{
      restPath=parsedEmuPath[2:length(parsedEmuPath)]
    }
    if(p1[['pattern']]){
      dirs=list.dirs(currentPath,full.names=FALSE,recursive=FALSE)
      bll=list()
      for(dir in dirs){
        dirPatt=gsub('*','.*',p1[['dir']],fixed=TRUE)
        if(grepl(dirPatt,dir)){
          newPath=file.path(currentPath,dir)
          bl=build_legacyBundleList(restPath,newPath,fileSuffixPattern,bundleList=bundleList)
          # prepend dir to list
          bl=lapply(bl,function(x,s) return(c(s,x)),dir)
          bll=c(bll,bl)
        }
        
      }
      return(bll)
    }else{
      if(is.null(currentPath)){
        newPath=p1[['dir']]
      }else{
        newPath=file.path(currentPath,p1[['dir']])
      }
      
      bll=build_legacyBundleList(restPath,newPath,fileSuffixPattern,bundleList=bundleList)
      return(bll)
    }
  }
  
}

convert_legacyBundleId<-function(legacybundleID){
  # takes character vector of legacy globpattern dirs and bundle name
  # and converts to session and bundle
  # examples: 
  # legacybundleID: "BLOCK10","SES1000","foo_42" -> "BLOCK10_SES1000","foo_42"
  # legacybundleID: "SES1000","foo_42" -> "SES1000","foo_42"
  # legacybundleID: "foo_42" -> "0000","foo_42"
  
  legacybundleIDLen=length(legacybundleID)
  globPatternCount=legacybundleIDLen-1
  
  if(globPatternCount>0){
    # collapse globpattern matches to one session ID
    s=paste(legacybundleID[1:globPatternCount],collapse='_')
    
  }else{
    # no glob patterns, put all bundles to dummy session
    s='0000'
    
  }
  return(c(s,legacybundleID[legacybundleIDLen]))
}

get_legacyEmuBundles=function(basePath,pathPattern,primaryFileSuffixPattern=NULL){
  if(is_relativeFilePath(pathPattern)){
    absPathPattern=file.path(basePath,pathPattern)
  }else{
    absPathPattern=pathPattern
  }
  emuParsedPathPattern=parse_emuTrackPath(absPathPattern)
  bl=build_legacyBundleList(emuParsedPathPattern[['dirs']],fileSuffixPattern=primaryFileSuffixPattern)
  return(bl)
}


convert_emuTrackPath <- function(absEmuTrackPath){
  # Emu track path may have asterisks for pattern matching e.g.
  # wav E:/KielCorpusRead/*/*/*/*
  # to use this pattern with R regex we have to convert to regular expression
  epSpl=strsplit(absEmuTrackPath,'/')[[1]]
  pp=c()
  cDir=NULL
  lastIsAsterisk=FALSE
  for(epDir in epSpl){
    #if(epDir!='*'){
    if(!grepl('*',epDir)){
      # Bug this condition block is never reached!!
      # The regex substitution is done is list.trackdirs
      lastIsAsterisk=FALSE

      if(epDir!=''){
        if(is.null(cDir)){
          cDir=epDir
        }else{
          sep='/'
          if(length(gl <- grep('/$', cDir))){
            sep=''
          }
          cDir=paste(cDir,epDir,sep=sep)
        }
      }else{
        cDir='/'
      }
    }else{
      lastIsAsterisk=TRUE
      
      pp=c(pp,cDir,epDir)
      cDir=NULL
    }
  }
  if(!lastIsAsterisk){
    pp=c(pp,cDir)
  }
  return(pp)
}

parse_emuTrackPath <- function(absEmuTrackPath){
  # Emu track path may have asterisks for pattern matching e.g.
  # wav E:/KielCorpusRead/*/*/*/*
  # to use this pattern with R regex we have to convert to regular expression
  epSpl=strsplit(absEmuTrackPath,'/+')[[1]]
  topo=list()
  wildcardDirLevelCount=0
  pp=list()
  cDir=NULL
  lastIsAsterisk=FALSE
  for(epDir in epSpl){
    if(!grepl('[*]',epDir)){
      lastIsAsterisk=FALSE
      if(epDir!=''){
        if(is.null(cDir)){
          cDir=epDir
        }else{
          sep='/'
          if(length(gl <- grep('/$', cDir))){
            sep=''
          }
          cDir=paste(cDir,epDir,sep=sep)
        }
      }else{
        cDir='/'
      }
    }else{
      lastIsAsterisk=TRUE
      
      if(!is.null(cDir)){
        pp[[length(pp)+1]]=list(dir=cDir,pattern=FALSE)
      } 
      cDir=NULL
      
      pp[[length(pp)+1]]=list(dir=epDir,pattern=TRUE)
      wildcardDirLevelCount=wildcardDirLevelCount+1
      
      
    }
  }
  if(!lastIsAsterisk){
    pp[[length(pp)+1]]=list(dir=cDir,pattern=FALSE)
  }
  topo[['dirs']]=pp
  topo[['patternCount']]=wildcardDirLevelCount
  
  return(topo)
}

## List directories to search for track files
## @description emu track pathes syntax allows asterisk wildcard pattern for path. This function goes through the directory hierarchy and returns list with all directories matching the pattern 
## @param emuPath Emu path specification (may contain asterisk wildcards) type character
## @param parsedEmuPathPattern character vector containing the parsed segments of the path. Each segment is a dierctory or a wildacrd asterisk.
## @return character vector of absolute path directories 
## @keywords emuDB bundle Emu
## 
list_trackdirs<-function(emuPath=NULL,parsedEmuPathPattern=NULL){
  if(is.null(parsedEmuPathPattern)){
    if(is.null(emuPath)){
      stop("At least one of the parameters emuPath or parsedEmuPathPattern is required.")
    }
    parsedEmuPathPattern=convert_emuTrackPath(emuPath)
  }
  cDir=NULL
  dirLevels=length(parsedEmuPathPattern)
  res=c()
  for(i in 1:dirLevels){
    ettp=parsedEmuPathPattern[i]
    lastLevel=(i==dirLevels)
    if(!grepl('[*]',ettp)){
      if(is.null(cDir)){
        cDir=ettp
      }else{
        cDir=file.path(cDir,ettp)
      }
      if(lastLevel){
        return(cDir)
      }
    }else{
      
      dirs=list.dirs(cDir,recursive=FALSE)
      for(dir in dirs){
        dirPatt=gsub('*','.*',ettp,fixed=TRUE)
        if(grepl(dirPatt,dir)){
          newPattern=c(dir)
          if(!lastLevel){
            for(j in (i+1L):dirLevels){
              newPattern=c(newPattern,parsedEmuPathPattern[j])
            }
          }
          wcRes=list_trackdirs(parsedEmuPathPattern=newPattern)
          res=c(res,wcRes)
        }
      }
      return(res)
    }
  }
  return(cDir)
}

get_legacyFilePath=function(basePath,emuPath,legacybundleID,fileExtension){
  if(is_relativeFilePath(emuPath)){
    absPathPattern=file.path(basePath,emuPath)
  }else{
    absPathPattern=emuPath
  }
  pp=parse_emuTrackPath(absEmuTrackPath = absPathPattern)
  path=NULL
  bIdIdx=1
  for(pdl in pp[['dirs']]){
    
    if(pdl[['pattern']]){
      # substitute
      dir=legacybundleID[bIdIdx]
      bIdIdx=bIdIdx+1
    }else{
      dir=pdl[['dir']]
    }
    if(is.null(path)){
      path=dir
    }else{
      path=file.path(path,dir)
    }
  }
  filename=paste0(legacybundleID[length(legacybundleID)],'.',fileExtension)
  return(file.path(path,filename))
}

## @import stringr wrassp
load_annotationForLegacyBundle=function(schema,legacyBundleID,basePath=NULL,encoding=NULL){
  
  newBundleId=convert_legacyBundleId(legacyBundleID)
  bundleName=newBundleId[2]
  sessionName=newBundleId[1]
  # determine samplerate
  # fallback is primary file
  sampleRateReferenceFile=NULL
  sampleTrackFile=NULL
  if(!is.null(schema[['mediafileBasePathPattern']]) && ! is.null(schema[['mediafileExtension']])){
    # use samples track to determine sample rate
    
    ## resolve wildcards
    #sampleRateReferenceFile=find.file.in.emu.path.pattern(emuPathPattern=schema[['mediafileBasePathPattern']],fileName=sampleTrackFile,basePath)
    sampleRateReferenceFile=get_legacyFilePath(basePath,emuPath=schema[['mediafileBasePathPattern']],legacyBundleID,fileExtension=schema[['mediafileExtension']])
  }
  if(is.null(sampleRateReferenceFile)){
    stop("Could not determine media sample rate of bundle ID ",paste(legacyBundleID,collapse='_'),"\n")
  }else{
    # TODO ASSP does not return good error messages if an IO error (not exist, permission dnied ,etc...) occurs
    # TODO test file access first
    pfAssp=read.AsspDataObj(sampleRateReferenceFile,0,4000)
    sampleRate=attr(pfAssp,'sampleRate')
  }
  
  # create signal paths list
  signalpaths=list() 
  for(tr in schema[['tracks']]){
    sFile=get_legacyFilePath(basePath=basePath,emuPath=tr[['basePath']],legacyBundleID,fileExtension=tr[['fileExtension']])
    if(!is.null(sFile)){
      signalpaths[[length(signalpaths)+1L]]=sFile
    }
  }
  idCnt=0
  levels=list()
  links=list()
  #  ESPS label files first
  for(ad in schema[['annotationDescriptors']]){
    extension=ad[['extension']]
    annoBasePath=NULL
    if(is.null(ad[['basePath']])){
      # TODO use same as primary track
    }else{
      annoBasePath=ad[['basePath']]
      # Emu: assume that files reside in this directory (no recursive search)
      annoPath=get_legacyFilePath(basePath=basePath,emuPath=ad[['basePath']],legacyBundleID,fileExtension=extension)
      if(!is.null(annoPath)){
        if(extension!='hlb'){
          # parse lab file
          if(file.exists(annoPath)){
            labTier=parse_espsLabelFile(labFilePath=annoPath,tierName=ad[['name']],tierType=ad[['type']],encoding=encoding,sampleRate=sampleRate,idCnt=idCnt)
            if(!is.null(labTier)){
              levels[[labTier[['name']]]] <- labTier
              labTierItemCnt=length(labTier[['items']])
              idCnt=idCnt+labTierItemCnt
            }
          }else{
            # warning ??
          }
        }
      }
    }
    
  }
  # now hlb file
  for(ad in schema[['annotationDescriptors']]){
    extension=ad[['extension']]

    if(is.null(ad[['basePath']])){
      # TODO use same as primary track
    }else{
      annoBasePathEmu=ad[['basePath']]
      # resolve wildcards
      if(extension=='hlb'){
        hlbFilePath=get_legacyFilePath(basePath=basePath,emuPath=annoBasePathEmu,legacyBundleID,fileExtension=extension)
        if(file.exists(hlbFilePath)){
          hlbParseResult=parse_hlbFile(hlbFilePath=annoPath,levelDefinitions=schema[['levelDefinitions']],levels=levels,encoding=encoding);
          hlbTiers=hlbParseResult[['hlbTiers']]
          links=hlbParseResult[['links']]
          # sort levels
          lIdx=0
          sortedLevels=list()
          
          for( ld in schema[['levelDefinitions']]){
            lIdx=lIdx+1L
            for(hlbTier in hlbTiers){
              if(ld[['name']]==hlbTier[['name']]){
                sortedLevels[[hlbTier[['name']]]] <- hlbTier
                break;
              }
            }
          }
          levels=sortedLevels
        }else{
          #cat("Warning: HLB file: ",hlbFilePath," does not exist!\n")
        }
      }
    }
    
  }
  
  bundleSampleRate=NULL
  for(l in levels){
    lvlSr=l[['sampleRate']]
    if(is.null(bundleSampleRate)){
      if(!is.null(lvlSr)){
        bundleSampleRate=lvlSr
      }
      
    }else{
      if(!is.null(lvlSr) && lvlSr!=bundleSampleRate){
        cat("WARNING: Levels have different sample rates!\n")
      }
    }
  }
  # set sample rate even if no annotation levels exist
  if(is.null(bundleSampleRate)){
    bundleSampleRate=sampleRate
  }
  
  
  sampleTrackFile=paste0(bundleName,'.',schema[['mediafileExtension']]) 
  annotates=paste0(sampleTrackFile)
  bundle=list(name=bundleName,sessionName=sessionName,legacyBundleID=legacyBundleID,annotates=annotates,sampleRate=bundleSampleRate,levels=levels,signalpaths=signalpaths,mediaFilePath=sampleRateReferenceFile,links=links)
  return(bundle)
}


build_hashedLinkDefs<-function(linkDefinitions){
  
  # build link definitions hashed by super level name
  linkDefsHashed=list()
  for(ld in linkDefinitions){
    supLvlNm=ld[['superlevelName']]
    linkDefsHashed[[supLvlNm]]=c(linkDefsHashed[[supLvlNm]],ld[['sublevelName']])
  }
  return(linkDefsHashed)
}

remove_redundantBundleLinks<-function(linkDefsHashed,bundle){
  lvls=bundle[['levels']]
  itemsHashed=list()
  for(lvl in lvls){
    for(it in lvl[['items']]){
      itemsHashed[[it[['id']]+1]]=lvl[['name']]
    }
  }
  legacyLinks=bundle[['links']]
  
  # new link list without redundant links
  links=list()
  for(legLk in legacyLinks){
    fromLvl=itemsHashed[[legLk[['fromID']]+1]]
    toLvl=itemsHashed[[legLk[['toID']]+1]]
    
    subLvls=linkDefsHashed[[fromLvl]]
    for(subLvl in subLvls){
      if(subLvl==toLvl){
        links[[length(links)+1]]=legLk
      }
    }
  }
  # overwrite
  bundle[['links']]=links
  return(bundle)
  
}


##' @title Convert legacy EMU database to the emuDB format
##' @description Converts an existing legacy EMU database to emuDB database structure. 
##' Copies or rewrites signal files and converts the database configuration and annotation data.
##' The legacy database must be addressed by its template file.
##' @details The database will be converted if the legacy database template file \code{emuTplPath} could 
##' be found and successfully loaded and parsed. The legacy template file usually has the extension '.tpl'. 
##' The UUID of the new emuDB will be randomly generated by default. If \code{targetDir} does not exist, 
##' the directory and its parents will be created. A new directory with the name of the database and the 
##' suffix '_emuDB' will be created in the \code{targetDir}. If the new database directory exists 
##' already, the function stops with an error. The template file is converted to a JSON file.
##' 
##' Some of the flags of the legacy EMU template files are ignored (lines with this syntax: "set [flagName] [flagValue]", 
##' known ignored flag names are: 'LabelTracks', 'SpectrogramWhiteLevel', 'HierarchyViewLevels', 'SignalViewLevels'). 
##' Legacy EMU utterances are reorganized to sessions and bundles. The naming of the sessions depends on the wildcard 
##' path pattern of the primary track: If the path contains no wildcard, only one session with the name '0000' will be created. 
##' If the path contains one wildcard path element, the names of the directories matching the pattern will be used as session names. 
##' If the path contains more than one wildcard path element, the session name is the concatenation of directory names 
##' separated by an underscore character.
##' 
##' Media files (usually WAV files) are copied, SSFF track files are rewritten using the ASSP library of package
##' \code{wrassp} by default (see option \code{rewriteSSFFTracks} below, see also \link[wrassp]{read.AsspDataObj} 
##' \link[wrassp]{write.AsspDataObj}). Annotations in EMU hierarchy (.hlb) files and ESPS label files are 
##' converted to one JSON file per bundle (utterance). Only those files get copied, which match the scheme 
##' of the template file. Additional files in the legacy database directories are ignored. The legacy EMU 
##' database will not be modified. For more information on the structural elements of an emuDB see \code{vignette{emuDB}}.
##' 
##'
##' \code{options} is a list of key value pairs:
##' @param emuTplPath EMU template file path
##' @param targetDir target directory
##' @param dbUUID optional UUID of emuDB, will be generated by default
##' @param ... currently available additional options:
##' \itemize{
##' \item{\code{rewriteSSFFTracks}: if \code{TRUE}, rewrite SSFF tracks instead of copying the file to get rid of big endian encoded SSFF files (SPARC), default: \code{TRUE}}
##' \item{\code{ignoreMissingSSFFTrackFiles}: if \code{TRUE}, missing SSFF track files are ignored, if \code{FALSE} an error will be generated, default: \code{TRUE}}
##' \item{\code{sourceFileTextEncoding}: encoding of legacy database text files (template, label and hlb files), possible values: NULL, "latin1", "UTF-8" "bytes" or "unknown" :default \code{NULL} (uses encoding of operating system platform)}
##' \item{\code{symbolicLinkSignalFiles}: if \code{TRUE}, signal files are symbolic linked instead of copied. Implies: \code{rewriteSSFFTracks=FALSE}, Default: \code{FALSE}}
##' }
##' @param verbose be verbose, default: \code{TRUE}
##' @seealso \code{\link{load_emuDB}}
##' @export
##' @name convert_legacyEmuDB
##' @keywords emuDB database schema Emu
##' @examples
##' \dontrun{
##' ## Convert legacy EMU database specified by EMU 
##' ## template file /mydata/EMU_legacy/ae/ae.tpl to directory /mydata/EMU/
##' ## and load it afterwards
##'
##' convert_legacyEmuDB("/mydata/EMU_legacy/ae/ae.tpl","/mydata/EMU/")
##' ae=load_emuDB("/mydata/EMU/ae_emuDB")
##'
##' ## Convert database "ae" and do not rewrite SSFF tracks 
##' 
##' convert_legacyEmuDB("/mydata/EMU_legacy/ae/ae.tpl",
##' "/mydata/EMU/",
##' options=list(rewriteSSFFTracks=FALSE))
##' 
##' ## Convert legacy database "ae" from emuR demo data and load converted emuDB
##' 
##' create_emuRdemoData()
##' demoTplPath=file.path(tempdir(),"emuR_demoData/legacy_ae/ae.tpl")
##' targetDir=file.path(tempdir(),"converted_to_emuR")
##' convert_legacyEmuDB(demoTplPath,targetDir)
##' dbHandle=load_emuDB(file.path(targetDir,"ae_emuDB"))
##' 
##' }
##' 
convert_legacyEmuDB <- function(emuTplPath,targetDir,dbUUID=UUIDgenerate(),verbose=TRUE, ...){
  # get ... options
  options = list(...)
  
  progress=0
  # default options
  # ignore missing SSFF track files
  # rewrite SSFF track files
  # encoding : platform 
  mergedOptions=list(sourceFileTextEncoding=NULL,ignoreMissingSSFFTrackFiles=TRUE,rewriteSSFFTracks=TRUE,symbolicLinkSignalFiles=FALSE)
  if(!is.null(options)){
    for(opt in names(options)){
      mergedOptions[[opt]]=options[[opt]]
    }
  }
  
  # pre check target dir
  if(file.exists(targetDir)){
    tdInfo=file.info(targetDir)
    if(!tdInfo[['isdir']]){
      stop(targetDir," exists and is not a directory.")
    }
  }
  legacyBasePath=dirname(emuTplPath)
  # load database schema and metadata to get db name
  dbConfig=load_dbConfigFromEmuTemplate(emuTplPath,dbUUID=dbUUID,encoding=mergedOptions[['sourceFileTextEncoding']])
  # database dir
  pp=file.path(targetDir,paste0(dbConfig[['name']],emuDB.suffix))
  
  # check existence of database dir
  if(file.exists(pp)){
    stop("Database storage dir ",pp," already exists.")
  }
  
  
  progress=progress+1L
  
  tplBaseDir=NULL
  tplBaseDir=dirname(emuTplPath)
  
  # create database dir in targetdir
  dir.create(pp,recursive = TRUE)
  
  # get UUID
  dbUUID=dbConfig[['UUID']]
  
  # get name
  dbName=dbConfig[['name']]
  
  # set user editable
  dbConfig[['EMUwebAppConfig']][['activeButtons']]=list(saveBundle=TRUE, showHierarchy=TRUE)

  # add handle for in memory DB
  dbHandle = emuDBhandle(dbName, pp, dbUUID, connectionPath = ":memory:")
  
  # filter transient properties
  # the properties listed are not persisted to JSON files
  
#   (filter list of old persitence code)  
#   emuR.persist.filters.DBconfig=list()
#   emuR.persist.filters.DBconfig[[1]]=c('annotationDescriptors')
#   emuR.persist.filters.DBconfig[[2]]=c('tracks')
#   emuR.persist.filters.DBconfig[[3]]=c('flags')
#   emuR.persist.filters.DBconfig[[4]]=c('ssffTrackDefinitions','basePath')
#   emuR.persist.filters.DBconfig[[5]]=c('mediafileBasePathPattern')
#   emuR.persist.filters.DBconfig[[6]]=c('maxNumberOfLabels')
#   emuR.persist.filters.DBconfig[[7]]=c('itemColNames')
#   emuR.persist.filters.DBconfig[[8]]=c('basePath')
#   emuR.persist.filters.DBconfig[[9]]=c('DBconfigPath')
  
  dbConfigPersist=dbConfig
  for(pNm in names(dbConfig)){
    if(pNm %in% c('annotationDescriptors','tracks','flags','mediafileBasePathPattern','maxNumberOfLabels','itemColNames','basePath','DBconfigPath')){
      dbConfigPersist[[pNm]]=NULL
    }
  }
  ssffTrCnt=length(dbConfig[['ssffTrackDefinitions']])
  if(ssffTrCnt>0){
    for(i in 1:ssffTrCnt){
      dbConfigPersist[['ssffTrackDefinitions']][[i]][['basePath']]=NULL
    }
  }
  
  # store db schema file
  store_DBconfig(dbHandle, dbConfigPersist)
  progress=progress+1L
  
  # load primary track file list first
  # and find samples track to get sample rate
  primaryFileList = NULL
  
  primaryBasePath=NULL
  primaryFileExtension=NULL
  
  primaryFileSuffixPattern=NULL
  # find primary and sample track paths
  for(tr in dbConfig[['tracks']]){
    if(tr[['fileExtension']]==dbConfig[['flags']][['PrimaryExtension']]){
      primaryFileExtension=tr[['fileExtension']]
      primaryBasePath=tr[['basePath']]
    }
  }
  
  if(is.null(primaryFileExtension)){
    for(ad in dbConfig[['annotationDescriptors']]){
      if(ad[['extension']]==dbConfig[['flags']][['PrimaryExtension']]){
        primaryFileExtension=ad[['extension']]
        primaryBasePath=ad[['basePath']]
        break
      }
    }
    
  }
  primaryFileSuffixPattern=paste0('[.]',primaryFileExtension,'$')
  legacyBundleIDsList=get_legacyEmuBundles(legacyBasePath,primaryBasePath,primaryFileSuffixPattern)
  
  
  bundlesCount=length(legacyBundleIDsList)

  us=1:bundlesCount
  if(verbose){
    cat("INFO: Converting legacy EMU database containing",bundlesCount,"bundles...\n")
    pb=utils::txtProgressBar(min=0,max=bundlesCount+2,initial=progress,style=3)
    
    utils::setTxtProgressBar(pb,progress)
  }
  linkDefsHashed=build_hashedLinkDefs(dbConfig[['linkDefinitions']])
  for(ui in us){
    legacyBundleID=legacyBundleIDsList[[ui]]
    newBundleId=convert_legacyBundleId(legacyBundleID)
    sessionName=newBundleId[1]
    bundleName=newBundleId[2]
    sDir=paste0(sessionName,session.suffix)
    sfp=file.path(pp,sDir)
    if(!file.exists(sfp)){
      dir.create(sfp)
    }
    ptrFilePath=get_legacyFilePath(legacyBasePath,primaryBasePath,legacyBundleID,primaryFileExtension)

    ptrFileBasename=basename(ptrFilePath)

    cutLen=str_length(primaryFileExtension)+1L
    cutPos=str_length(ptrFileBasename)-cutLen

    bundle=load_annotationForLegacyBundle(dbConfig,legacyBundleID,legacyBasePath,encoding=mergedOptions[['sourceFileTextEncoding']])
    bundle=remove_redundantBundleLinks(linkDefsHashed,bundle)

    bDir=paste0(bundle[['name']],bundle.dir.suffix)
    bfp=file.path(sfp,bDir)
    dir.create(bfp)
    # create new list that only contains relevant infos
    bp = list(name = bundle$name, annotates = bundle$annotates, 
              sampleRate = bundle$sampleRate, levels = bundle$levels, links = bundle$links)
    # remove sample rate entries
    for(i in 1:length(bp$levels)){
      bp$levels[[i]]$sampleRate = NULL
    }
    # remove level names
    names(bp$levels) = NULL
    
    # metadata (annotations)
    ban=str_c(bundle[['name']],bundle.annotation.suffix,'.json')
    baJSONPath=file.path(bfp,ban)
    pbpJSON=jsonlite::toJSON(bp,auto_unbox=TRUE,force=TRUE,pretty=TRUE)
    writeLines(pbpJSON,baJSONPath)
    
    
    for(sf in bundle[['signalpaths']]){
      #cat("Signalpath: ",sf,"\n")
      bn=basename(sf)
      nsfp=file.path(bfp,bn)
      # check if SSFF type
      isSSFFFile=FALSE
      for(ssffTrDef in dbConfig[['ssffTrackDefinitions']]){
        ssffTrFileExt=ssffTrDef[['fileExtension']]
        fileExtPatt=paste0('[.]',ssffTrFileExt,'$')
        if(length(grep(fileExtPatt,sf))==1){
          isSSFFFile=TRUE
          break
        }
      }
      if(file.exists(sf)){
        if(mergedOptions[['symbolicLinkSignalFiles']]){
          file.symlink(from=sf,to=nsfp)
        }else if(mergedOptions[['rewriteSSFFTracks']] && isSSFFFile){
          # is SSFF track
          # read/write instead of copy to get rid of big endian encoded SSFF files (SPARC)
          pfAssp=read.AsspDataObj(sf)
          write.AsspDataObj(pfAssp,nsfp)
        }else{
          # media file (likely a wav file)
          file.copy(from=sf,to=nsfp)
        }
      }else{
        if(!mergedOptions[['ignoreMissingSSFFTrackFiles']]){
          stop("SSFF track file :'",sf,"' does not exist!")
        }
      }
    }
    bundle[['levels']]=NULL
    bundle[['links']]=NULL
    
    bName=bundle[['name']]
    
    progress=progress+1L
    if(verbose){
      utils::setTxtProgressBar(pb,progress)
    }
  }
  if(verbose){
    utils::setTxtProgressBar(pb,progress)
    cat("\n")
  }
}
