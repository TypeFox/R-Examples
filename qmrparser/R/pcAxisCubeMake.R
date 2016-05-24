#do not edit, edit noweb/qmrparser.nw
  pcAxisCubeMake <- function(cstream) {
##
    stringsAsFactorsOld <- options("stringsAsFactors")
    on.exit(options(stringsAsFactors=stringsAsFactorsOld$stringsAsFactors
))
    options(stringsAsFactors=FALSE)
##
##      Creates a vector with values associated with a keyword, taking into account data type
    avecMapValue <- function(e) {
      mapValue <- e[[4]]
      if (length(mapValue) == 0 ) return(c())
      switch(mapValue$type,
             stringstring=as.character(mapValue$value),
             liststring  =as.character(mapValue$value),
             symbol      =as.character(mapValue$value),
             list        =as.character(mapValue$value),
             tlist       =c(as.character(unlist(mapValue$value[[1]]$value)),as.character(mapValue$value[[2]]$type),as.character(unlist(mapValue$value[[2]]$value))),
             as.character(mapValue$value))
    }
##
##      Flattens data associated with a keyword
##      "rownum"  "keyword"  "variable"  "variableValue"  "mapValueType"   "mapValueLength" "mapValue"

    keywordsAplana <- function(e)
      list(   keyword = e[[1]][[1]],
           language= if(! is.null(e[[2]]) && length(e[[2]]) > 0 )  e[[2]][[1]] else list(as.character("")),
           arity   =length(e[[3]]),
           args    = e[[3]],
           type    = e[[4]]$type,
           length  = length(e[[4]]$value),
           value   = avecMapValue(e)       )
    ##
    ##
    productoCartersiano <- function(a,b) unlist(lapply(a,function(e1) lapply(b,function(e2) c(list(),e1,e2))),recursive=FALSE,use.names = FALSE)
    ##  
    combinaValores <- function(listVar) { 
      listVar <- Filter(function(x) length(x) > 0, listVar)
      if( length(listVar) <= 1 ) as.list(listVar[[1]]) 
      else { 
        expa <- listVar[[1]] 
        for( i in 2:length(listVar) ) expa <- productoCartersiano(expa,listVar[[i]])
        return(expa) 
      }
    }
    proyect <- function(listlist,n) lapply(listlist,function(e) e[[n]])
    ##
    ##
    keywords     <- lapply(cstream[[2]],function(x) keywordsAplana(x))

    selectClauses <- function(data,filter,fields)
      if( missing(filter) ) { if ( missing(fields) )               data  else lapply(              data ,function(e) e[fields]) } else 
    { if ( missing(fields) ) Filter(filter,data) else lapply(Filter(filter,data),function(e) e[fields]) }

    getColumn     <- function(data,field)             as.vector(unlist(lapply(data,function(e) e[[field]])))
    ##
    STUB          <- selectClauses(keywords, function(e) e$language=="" && e$arity == 0 && e$keyword == "STUB"   )[[1]]
    HEADING       <- selectClauses(keywords, function(e) e$language=="" && e$arity == 0 && e$keyword == "HEADING")[[1]]
    DATA          <- selectClauses(keywords, function(e) e$language=="" && e$arity == 0 && e$keyword == "DATA"   )[[1]]
    TIMEVAL       <- selectClauses(keywords, function(e) e$language=="" && e$arity == 1 && e$keyword == "TIMEVAL")
    KEYS          <- selectClauses(keywords, function(e) e$language=="" && e$arity == 1 && e$keyword == "KEYS"   )
    ##
    frequency      <- if ( length(TIMEVAL) > 0 ) TIMEVAL[[1]]$value[[1]] else as.character(NA) 
    keysFlag       <-      length(KEYS)    > 0 
    ##
    ## pxCube:
    ## (headingLength, StubLength, frequency)
    pxCube <-data.frame(headingLength=HEADING$length,StubLength=STUB$length,frequency=frequency)
    ##
    ## pxCubeAttrN
    aritys       <- unique(as.vector(lapply(keywords,function(e) e$arity)))
    pxCubeAttrN  <- list() 
    for ( a in aritys ) {
      ##print(a)
      ks          <- Filter(function(e) e$arity == a && e$keyword != "DATA" , keywords )         
      ksap        <- lapply(ks, function(e) c(e$keyword,e$language,e$args,e$length,value=paste('"',e$value,'"',sep="",collapse=" ")[[1]]))
      ksdf        <- data.frame(lapply(do.call(function(...) rbind.data.frame(..., deparse.level=0),ksap),unlist))
      names(ksdf) <- c("keyword","language",if(a>0)paste("arg",1:a,sep="") else NULL,"length","value")
      pxCubeAttrN[[paste("A",as.character(a),sep="")]] <-  ksdf         
    }
    ##
    ##
    ##
    ## CODES or VALUES registries
    codesvalues   <- selectClauses(keywords,
                                   function(e) e$language=="" && e$arity == 1 && e$keyword %in% c("CODES","VALUES"),
                                   c("keyword","args","length","value"))
    ##  Value number, by variables
    codesValuesLength       <- aggregate(getColumn(codesvalues,"length"),list(getColumn(codesvalues,"args")),max)
    rownames(codesValuesLength) <- codesValuesLength$Group.1
    colnames(codesValuesLength) <- c("variable","valueLength")
    ## Whether use "codes" or "values" in the numerical value key/pk
    codesOValues      <- aggregate(getColumn(codesvalues,"keyword"),list(getColumn(codesvalues,"args")),function(x) min(as.character(x)))
    rownames(codesOValues     ) <- codesOValues$Group.1
    colnames(codesOValues     ) <- c("variable","keyword")
    
    ## Information union: Number of values, by variables and whether use "codes" or "values" in the numerical value key/pk
    variableMeta           <- merge(codesOValues,codesValuesLength)
    rownames(variableMeta) <- codesOValues$variable
    ##

    valueselimination <-  selectClauses(keywords,function(e) e$language=="" && e$arity == 1 && e$keyword == "ELIMINATION")
    
    
    ## pxCubeVariable: 
    ## (variableName, headingOrStud, codesYesNo, valuesYesNo, variableOrder, valueLength)
    variableName  <- unlist(c(HEADING$value,STUB$value)) 
    
    headingOrStud <- c(rep("HEADING",HEADING$length),rep("STUB"   ,STUB$length)) 
    
    codesYesNo    <- sapply(variableName,function(v) length(selectClauses(codesvalues,function(e) e$keyword == "CODES"  && e$args[[1]] == v)) >= 1 )
    
    valuesYesNo   <- sapply(variableName,function(v) length(selectClauses(codesvalues,function(e) e$keyword == "VALUES" && e$args[[1]] == v)) >= 1 )
    
    variableOrder <- c(1:HEADING$length,1:STUB$length)
    
    valueLength   <- unlist(variableMeta[variableName,"valueLength"])
    
    pxCubeVariable<- cbind.data.frame(variableName,headingOrStud, codesYesNo, valuesYesNo, variableOrder, valueLength)
    
    colnames(pxCubeVariable) <- c("variableName","headingOrStud","codesYesNo","valuesYesNo","variableOrder","valueLength")
    ##
    ## pxCubeVariableDomain:
    ## (variableName, code, value, valueOrder, eliminationYesNo)
    variableName     <- unlist(sapply(1:length(pxCubeVariable$variableName),function(i) rep(pxCubeVariable$variableName[i],pxCubeVariable$valueLength[i])))
    
    code             <- unlist(sapply(1:length(pxCubeVariable$variableName),function(i) if(pxCubeVariable$codesYesNo[i])                               selectClauses(codesvalues,function(e) e$keyword == "CODES"  && e$args[[1]] == pxCubeVariable$variableName[i],"value") else rep(as.character(NA),pxCubeVariable$valueLength[i])))
    
    value                <- unlist(sapply(1:length(pxCubeVariable$variableName),function(i) if(pxCubeVariable$valuesYesNo[i])                       selectClauses(codesvalues,function(e) e$keyword == "VALUES" && e$args[[1]] == pxCubeVariable$variableName[i],"value") else rep(as.character(NA),pxCubeVariable$valueLength[i])))
    
    valueOrder       <- unlist(sapply(1:length(pxCubeVariable$variableName),function(i) 1:pxCubeVariable$valueLength[i]))
    
    eliminationYesNo <- unlist(sapply(1:length(value),function(i) length(selectClauses(valueselimination,function(e) e$keyword == "ELIMINATION" && e$args[[1]] == variableName[i] && e$value == value[i])) > 0))

    if ( length(code) != length(value) ) stop("Error: CODE and VALUE non-consistent")
    
    pxCubeVariableDomain        <-  cbind.data.frame(variableName,code,value,valueOrder,eliminationYesNo)
        names(pxCubeVariableDomain) <- c("variableName","code","value","valueOrder","eliminationYesNo")
    ##
    ## pxCubeData:
    ## ({variableName}+, data)
    if ( ! keysFlag ) { 
      STUDValues           <- combinaValores(lapply(STUB$value   ,function(v) unlist(selectClauses(codesvalues,function(e) e$keyword == variableMeta[v,"keyword"] && e$args[[1]] == v,"value"))))
      
      HEADINGValues        <- combinaValores(lapply(HEADING$value,function(v) unlist(selectClauses(codesvalues,function(e) e$keyword == variableMeta[v,"keyword"] && e$args[[1]] == v,"value"))))
      
      pxCubeData           <- combinaValores(list(STUDValues,HEADINGValues))
      ##      
      ## data.frame construction, from the values combinations list
      e1         <- pxCubeData[[1]]
      pxCubeData <- do.call("cbind.data.frame",
                            lapply(1:length(e1),function(i) as.character(proyect(pxCubeData,i))))
      colnames(pxCubeData) <-  pxCubeVariable$variableName
      ##
      ##      
      if( length(DATA$value) < nrow(pxCubeData) ) stop("Error: variables and data inconsistency"," ",length(DATA$value) ," ", nrow(pxCubeData) )                
      
      if( length(DATA$value) > nrow(pxCubeData) ) warning("Warnings: variables and data inconsistency"," ",length(DATA$value) ," ", nrow(pxCubeData) )                
      
      pxCubeData$data      <- DATA$value[1:nrow(pxCubeData)]
      names(pxCubeData)    <- c(STUB$value,HEADING$value,"data")
      rownames(pxCubeData) <- 1:nrow(pxCubeData)

    } else {
      HEADINGValues        <- combinaValores(lapply(HEADING$value,function(v) unlist(selectClauses(codesvalues,function(e) e$keyword == variableMeta[v,"keyword"] && e$args[[1]] == v,"value"))))
      
      numFields              <- pxCube$StubLength + length(HEADINGValues)
      keysdata               <- DATA$value
      ## print(paste(pxCube$StubLength," - ",length(HEADINGValues)))
      ## print(paste( length(keysdata)," - ",numFields))                
      numreg                 <- length(keysdata)/numFields
      
      keys                   <- sapply(1:numreg,function(i) as.list(keysdata[((i-1)*numFields+1):((i-1)*numFields+pxCube$StubLength)]),simplify = FALSE)
      
      pxCubeData           <- combinaValores(list(keys,HEADINGValues))
      ##
      ## data.frame construction, from the values combinations list
      e1         <- pxCubeData[[1]]
      pxCubeData <- do.call("cbind.data.frame",
                            lapply(1:length(e1),function(i) as.character(proyect(pxCubeData,i))))
      colnames(pxCubeData) <-  pxCubeVariable$variableName
      ##
      ##      
      data                   <- as.vector(sapply(1:numreg,function(i) keysdata[((i-1)*numFields+1+pxCube$StubLength):((i-1)*numFields+pxCube$StubLength+length(HEADINGValues))],simplify = TRUE))

      if( length(data) != nrow(pxCubeData) ) stop("Error: variables and data inconsistency"," ",length(data), " ",nrow(pxCubeData))
      
      pxCubeData$data      <- as.numeric(data)
      names(pxCubeData)    <- c(STUB$value,HEADING$value,"data")
      rownames(pxCubeData) <- 1:nrow(pxCubeData)
      
    }
    ##
    return(list(pxCube=pxCube,pxCubeAttrN=pxCubeAttrN,pxCubeVariable=pxCubeVariable,pxCubeVariableDomain=pxCubeVariableDomain,pxCubeData=pxCubeData))
        
}
