requireNamespace("stringr", quietly = T)

## Parser for ESPS label files
## 
## @param labFilePath ESPS label file path
## @param tierName name of the tier
## @param sampleRate sample rate of corresponding signal file
## @param encoding text encoding (default: NULL -> R encoding "unknown")
## @return new tier containing parsed items
## @import stringr
## @keywords emuR ESPS lab Emu
## 
parse_espsLabelFile <- function(labFilePath=NULL,tierName,tierType=NULL,sampleRate,encoding=NULL,idCnt=0) {
  SIGNAL_KEY="signal"
  NUMBER_OF_FIELDS_KEY="nfields"
  SEPARATOR_KEY="separator"
  COLOR_KEY="color"
  DATA_SECTION_START_KEY='#'
  INTERVAL_FLAG_VALUE='H#'
  
  fileToRead=NULL
  inHeaderSection=TRUE
  intervalMode=FALSE
  firstDataLine=TRUE
  intervalStart=NULL
  itemList=list()

  if(is.null(labFilePath)){
      stop("Argument labFilepath or labCon must not be NULL\n")
  }else{
      fileToRead=labFilePath
  }  
  if(sampleRate <=0 ){
    stop("Samplerate must be greater than zero\n")
  }
  
    # read
    if(is.null(encoding)){
      lc = try(readLines(fileToRead))
    }else{
      lc = try(readLines(fileToRead,encoding=encoding))
    }
    if(class(lc) == "try-error") {
      stop("read.TextGrid: cannot read from file ", fileToRead)
    }
    
    for(line in lc){
      trimmedLine=stringr::str_trim(line)
      if(inHeaderSection){
      if(trimmedLine == DATA_SECTION_START_KEY){
        inHeaderSection=FALSE
      }else{
        kv=parse_lineToKeyValue(trimmedLine,'[[:space:]]')
        if(!is.null(kv) && kv[1]==NUMBER_OF_FIELDS_KEY){
          if(kv[2] != 1){
            stop("only files with one field supported")
          }
        }
        # ignore other headers
      }
    }else{
      lineTokensLst=str_split(trimmedLine,'[[:space:]]+',3)
      lineTokens=lineTokensLst[[1]]
      lineTokenCount=length(lineTokens);
      if(lineTokenCount>=2){
        
        timeStampStr=lineTokens[1]
        color=lineTokens[2]
        label=NULL
        if(lineTokenCount==3){
          label=lineTokens[3]
        }
        
        timeStamp=as(timeStampStr,"numeric")
        timeStampInSamples=timeStamp*sampleRate
        if(firstDataLine){
          
          if(is.null(tierType)){
            if(label==INTERVAL_FLAG_VALUE){
              tierType='SEGMENT'
            }else{
              tierType='EVENT'
            }
          }
        }
        
        if(firstDataLine & tierType=='SEGMENT'){
            samplePoint=floor(timeStampInSamples)
            intervalStartPoint=samplePoint
        }else{
          labelAttrs=list(list(name=tierName,value=label))
          #labelAttrs[[tierName]]=label
          
          if(tierType=='SEGMENT'){
            samplePoint=floor(timeStampInSamples)
            # duration calculation according to partitur format
            # the sum of all durations is not equal to the complete sample count
            duration=samplePoint-intervalStartPoint-1
            currItem=list(id=idCnt,sampleStart=intervalStartPoint,sampleDur=duration,labels=labelAttrs)
            idCnt=idCnt+1
            itemList[[length(itemList)+1]] <- currItem
            intervalStartPoint=samplePoint
            
          }else{
            samplePoint=round(timeStampInSamples)
            currItem=list(id=idCnt,samplePoint=samplePoint,labels=labelAttrs)
            idCnt=idCnt+1
            itemList[[length(itemList)+1]] <- currItem
          }
          
        }
        firstDataLine=FALSE
      }
    }
    }
  labTier=list(name=tierName,type=tierType,sampleRate=sampleRate,items=itemList);
  return(labTier)
}
