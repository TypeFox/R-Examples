requireNamespace("httpuv", quietly = T)
requireNamespace("base64enc", quietly = T)

.server_env = new.env()
assign("serverHandle", NULL, envir = .server_env)

getServerHandle <- function() {
  get("serverHandle", envir = .server_env)
}

setServerHandle <- function(sh) {
  assign("serverHandle", sh, envir = .server_env)
}

##' Serve EMU database to EMU-webApp
##' 
##' @description Serves emuDB media files, SSFF tracks and annotations for EMU-webApp browser GUI \url{http://ips-lmu.github.io/EMU-webApp/}
##' 
##' Instructions:
##' 
##' Start and connect:
##' 
##' \itemize{
##' \item Call this function to start the server.
##' \item Start a suitable HTML5 capable Web-Browser (Google Chrome, Firefox,...).
##' \item Navigate to the EMU-Webapp URL: \url{http://ips-lmu.github.io/EMU-webApp/}.
##' \item Press the 'Connect' button in the EMU-webApp and connect with default URL.
##' \item EMU-webApp loads the bundle list and the first bundles media file, SSFF tracks and annotations.
##' }
##' 
##' Disconnect and stop:
##' \itemize{
##' \item Disconnect and stop the server with the 'Clear' button of the webapp or the reload button of your browser.
##' \item The server can also be interrupted with Ctrl-C if something wents wrong.
##' \item To serve only a subset of sessions or bundles use the parameters \code{sessionPattern} and/or \code{bundlePattern}.
##' }
##' 
##' @details  Function opens a HTTP/websocket and waits in a loop for browser requests. Parameter host determines the IP address(es) of hosts allowed to connect to the server. By default the server only listens to localhost. If you want to allow connection from any host set the host parameter to \code{0.0.0.0}. Please note that this might be an safety issue! The \code{port} parameter determines the port the server listens on. The \code{host} and \code{port} parameters are intended only for expert users. When started the R console will be blocked. On successfull connection the server sends the session and bundle list of the database referenced by name by parameter \code{dbName} or by UUID parameter \code{dbUUID}.
##' The Web application requests bundle data for viewing or editing. If a bundle is modified with the EMU-webApp and the save button is pressed the server modifies the internal database and saves the changes to disk.
##' Communication between server and EMU webApp is defined by EMU-webApp-websocket-protocol version 0.0.2.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param sessionPattern A regular expression pattern matching session names to be served
##' @param bundlePattern A regular expression pattern matching bundle names to be served
##' @param host host IP to listen to (default: 127.0.0.1  (localhost))
##' @param port the port number to listen on (default: 17890)
##' @param autoOpenURL URL passed to \code{\link{browseURL}} function. If NULL or an empty string are passed in
##' \code{\link{browseURL}} will not be invoked.
##' @param debug TRUE to enable debugging (default: no debugging messages)
##' @param debugLevel integer higher values generate more detailed debug output
##' @return TRUE if the database was modified, FALSE otherwise
##' @export
##' @keywords emuDB EMU-webApp database websocket Emu
##' @examples
##' \dontrun{ 
##' ## Load EMU database 'myDb' and serve it to the EMU-webApp (opens default HTTP/websocket port 17890)
##' 
##' myDb = load_emuDB("/path/to/myDb")
##' serve(myDb)
##' }
##' 
serve <- function(emuDBhandle, sessionPattern='.*',bundlePattern='.*',host='127.0.0.1',port=17890, autoOpenURL = "http://ips-lmu.github.io/EMU-webApp/?autoConnect=true",  debug=FALSE,debugLevel=0){
  if(debug && debugLevel==0){
    debugLevel=2
  }
  modified=FALSE
  emuDBserverRunning=FALSE
  bundleCount=0
  DBconfig = load_DBconfig(emuDBhandle)

  allBundlesDf=list_bundles(emuDBhandle)
  bundlesDf=allBundlesDf
  if(!is.null(sessionPattern) && sessionPattern!='.*'){
    ssl=emuR_regexprl(sessionPattern,bundlesDf[['session']])
    bundlesDf=bundlesDf[ssl,]
  }
  if(!is.null(bundlePattern) && bundlePattern!='.*'){
    bsl=emuR_regexprl(bundlePattern,bundlesDf[['name']])
    bundlesDf=bundlesDf[bsl,]
  }
  
  
  httpRequest = function(req){
    # Only 
    # Rook conform answer  
    body = paste('<p>http protocol not supported, please use ws protocol.</p>')
    list(
      status = 501L,
      headers = list(
        'Content-Type' = 'text/html'
      ),
      body = body
    )
  }
  
  onHeaders<-function(req){
    # following httuv docs we shoul return NULL here to proceed but that terminates the R session!
    #return(NULL)
  }
  
  serverEstablished = function(ws){
    
    cat("emuR websocket service established\n")
    
    serverClosed = function(ws){
      
      cat("emuR websocket service closed\n")
      emuRserverRunning<<-FALSE
      
    }
    
    sendError = function(ws,errMsg,callbackID){
      status=list(type='ERROR',details=errMsg);
      response=list(callbackID=callbackID,status)
      responseJSON=jsonlite::toJSON(response,auto_unbox=TRUE,force=TRUE,pretty=TRUE)
      result=ws$send(responseJSON)
    } 
    
    serverReceive = function(isBinary,DATA){
      if(debugLevel >= 4 ){
        cat("onMessage() call, binary:",isBinary," data: ",DATA,"\n")
        
      }
      D = ""
      if(is.raw(DATA)) {
        D = rawToChar(DATA)
      }else{
        D = DATA
      }
      jr=jsonlite::fromJSON(D,simplifyVector = FALSE)
      if(debugLevel >= 2 ){
        cat("Received command from EMU-webApp: ",jr[['type']],"\n")
        if(debugLevel >= 3){
          jrNms=names(jr)
          
          for( jrNm in jrNms){
            value=jr[[jrNm]]
            cat("param: ",jrNm)
            if(class(value)=='character'){
              cat(": ",jr[[jrNm]])
            }
            cat("\n")
          }
        }
        
      }
      if(!is.null(jr$type)){
        if(debugLevel >= 2 ){
          cat("Received type from EMU-webApp: ",jr[['type']],"\n")
        }
        
      }
      if(jr$type == 'GETPROTOCOL'){
        
        protocolData=list(protocol='EMU-webApp-websocket-protocol',version='0.0.2')
        response=list(status=list(type='SUCCESS'),callbackID=jr$callbackID,data=protocolData)
        responseJSON=jsonlite::toJSON(response,auto_unbox=TRUE,force=TRUE,pretty=TRUE) 
        result=ws$send(responseJSON)
        if(debugLevel >= 2){
          cat("Sent protocol. \n")
        }
        
      }else if(jr$type == 'GETDOUSERMANAGEMENT'){
        # R server mode is single user mode 
        response=list(status=list(type='SUCCESS'),callbackID=jr$callbackID,data="NO")
        responseJSON=jsonlite::toJSON(response,auto_unbox=TRUE,force=TRUE,pretty=TRUE) 
        result=ws$send(responseJSON)
        if(debugLevel >= 2){
          cat("Sent user managment: no. \n")
        }
        
      }else if(jr$type == 'GETGLOBALDBCONFIG'){
        if(debugLevel >= 4){
          cat("Send config: ",as.character(DBconfig),"\n")
        }
        response=list(status=list(type='SUCCESS'),callbackID=jr$callbackID,data=DBconfig)
        responseJSON=jsonlite::toJSON(response,auto_unbox=TRUE,force=TRUE,pretty=TRUE) 
        result=ws$send(responseJSON)
        if(debugLevel >= 2){
          if(debugLevel >=4){
            cat(responseJSON,"\n")
          }
          cat("Sent config. \n")
        }
        #}
        
        
      }else if(jr$type == 'GETBUNDLELIST'){
        
        response=list(status=list(type='SUCCESS'),callbackID=jr$callbackID,dataType='uttList',data=bundlesDf)
        responseJSON=jsonlite::toJSON(response,auto_unbox=TRUE,force=TRUE,pretty=TRUE)
        
        if(debugLevel >= 5)cat(responseJSON,"\n")
        result=ws$send(responseJSON)
        if(debugLevel >= 2){
          cat("Sent utterance list with length: ",nrow(bundlesDf)," \n")
        }
        
      }else if(jr$type == 'GETBUNDLE'){
        
        bundleName=jr[['name']]
        bundleSess=jr[['session']]
        #cat("data:",jr[['data']],"\n")
        if(debugLevel>2){
          cat("Requested bundle:",bundleName,",session:",bundleSess,"\n")
        }
        err=NULL
        if(debugLevel>3){
          cat("Convert bundle to S3 format",bundleName,"\n")
        }
        # construct path to annotJSON
        annotFilePath = normalizePath(file.path(emuDBhandle$basePath, paste0(bundleSess, session.suffix), 
                                                paste0(bundleName, bundle.dir.suffix), 
                                                paste0(bundleName, bundle.annotation.suffix, '.json')))
        
        b = jsonlite::fromJSON(annotFilePath, simplifyVector = F)
        if(is.null(b)){
          # error
          err=simpleError(paste('Could not load bundle ',bundleName,' of session ',bundleSess))
        }
        if(is.null(err)){
          mediaFilePath=normalizePath(file.path(emuDBhandle$basePath, paste0(bundleSess, session.suffix), 
                                                paste0(bundleName, bundle.dir.suffix), 
                                                paste0(bundleName, ".", DBconfig$mediafileExtension)))
          if(debugLevel>4){
            cat("Mediafile: ",mediaFilePath," for ",b$name,"\n")
          }
          audioFile=tryCatch(file(mediaFilePath, "rb"),error=function(e) err<<-e)
          if(is.null(err)){
            audioFileData=readBin(audioFile, raw(), n=file.info(mediaFilePath)$size)
            if(inherits(audioFileData,'error')){
              err=audioFileData
            }else{
              audioBase64=base64enc::base64encode(audioFileData)
              mediaFile=list(encoding="BASE64",data=audioBase64)
              close(audioFile)
            }
          }
        }
        if(is.null(err)){   
          ssffTracksInUse=get_ssffTracksUsedByDBconfig(DBconfig)
          ssffTrackNmsInUse=c()
          for(ssffTrackInUse in ssffTracksInUse){
            ssffTrackNmsInUse=c(ssffTrackNmsInUse,ssffTrackInUse[['name']])
          }          
          if(debugLevel >= 4){
            
            cat(length(ssffTrackNmsInUse)," track definitions in use:\n")
            for(sfInU in ssffTrackNmsInUse){
              cat(sfInU," ")
            }
            cat("\n")
          }
          ssffFiles=list()
          # Hash (here: named character vector) with SSFF files extension as key and file path as value
          # avoids duplicates in ssff files list
          ssffFilesHash=character(0)
          for(ssffTr in DBconfig$ssffTrackDefinitions){
            if(ssffTr[['name']] %in% ssffTrackNmsInUse){
              fe=ssffTr[['fileExtension']]
              ssffFilesHash[fe]=normalizePath(file.path(emuDBhandle$basePath, paste0(bundleSess, session.suffix), 
                                                        paste0(bundleName, bundle.dir.suffix), 
                                                        paste0(bundleName, ".", fe)))
            }
          }
          # read SSFF track file data
          ssffFileExts=names(ssffFilesHash)
          for(ssffFileExt in ssffFileExts){
            ssffFilePath=ssffFilesHash[ssffFileExt]
            mf=tryCatch(file(ssffFilePath, "rb"),error=function(e) {err<<-e})
            if(is.null(err)){
              mfData=readBin(mf, raw(), n=file.info(ssffFilePath)$size)
              if(inherits(mfData,'error')){
                err=mfData
                break
              }
            }else{
              break
            }
            mfDataBase64=base64enc::base64encode(mfData)
            encoding="BASE64"
            ssffDatObj=list(encoding=encoding,data=mfDataBase64,fileExtension=ssffFileExt)
            ssffFiles[[length(ssffFiles)+1]]=ssffDatObj
            close(mf)
          }
          if(is.null(err)){
            data=list(mediaFile=mediaFile,ssffFiles=ssffFiles,annotation=b)
          }
        }
        
        if(is.null(err)){
          responseBundle=list(status=list(type='SUCCESS'),callbackID=jr$callbackID,responseContent='bundle',contentType='text/json',data=data)
        }else{
          errMsg=err[['message']]
          cat("Error: ",errMsg,"\n")
          responseBundle=list(status=list(type='ERROR',message=errMsg),callbackID=jr[['callbackID']],responseContent='status',contentType='text/json')
          
        }
        responseBundleJSON=jsonlite::toJSON(responseBundle,auto_unbox=TRUE,force=TRUE,pretty=TRUE)
        result=ws$send(responseBundleJSON)
        if(is.null(err) & debugLevel >= 2){
          
          if(debugLevel >=8){
            cat(responseBundleJSON,"\n")
          }
          cat("Sent bundle containing",length(ssffFiles),"SSFF files\n")
        }
        # reset error
        err=NULL
        
      }else if(jr[['type']] == 'SAVEBUNDLE'){
        jrData=jr[['data']]
        jrAnnotation=jrData[['annotation']]
        bundleSession=jrData[['session']]
        bundleName=jrData[['annotation']][['name']]
        if(debugLevel>3){
          cat("Save bundle ",bundleName," from session ",bundleSession,"\n");
        }
        err=NULL
        
        ssffFiles=jr[['data']][['ssffFiles']]
        oldBundleAnnotDFs = load_bundleAnnotDFsDBI(emuDBhandle, bundleSession, bundleName)
        
        # warnings as errors
        warnOptionSave=getOption('warn')
        options('warn'=2)
        responseBundle=NULL
        
        # check if cached version of bundle is available
        if(is.null(oldBundleAnnotDFs)){
          # error
          err=simpleError(paste('Could not load bundle ',bundleSession, bundleName))
        }else{
          for(ssffFile in ssffFiles){
            inCfg=FALSE
            sp = normalizePath(file.path(emuDBhandle$basePath, paste0(bundleSession, session.suffix), 
                                         paste0(bundleName, bundle.dir.suffix), 
                                         paste0(bundleName, ".", ssffFile$fileExtension)))
            if(is.null(sp)){
              errMsg=paste0("SSFF track definition for file extension '",ssffFile[['fileExtension']],"' not found!")
              err=simpleError(errMsg)
            }else{
              # store
              if(debugLevel>3){
                cat("Writing SSFF track to file: ",sp,"\n")
              }
              ssffTrackBin=base64enc::base64decode(ssffFile[['data']])
              ssffCon=tryCatch(file(sp,'wb'),error=function(e){err<<-e})
              if(is.null(err)){
                res=tryCatch(writeBin(ssffTrackBin,ssffCon))
                close(ssffCon)
                if(inherits(res,'error')){
                  err=res
                  break
                }
                modified<<-TRUE
              }
            }
          }
          bundleData=jr[['data']][['annotation']]

          # if we do not have an (error) response already
          if(is.null(err)){
            ##### emuDB ####
            # construct path to annotJSON and store
            annotFilePath = file.path(emuDBhandle$basePath, paste0(bundleSession, session.suffix), 
                                      paste0(bundleName, bundle.dir.suffix), 
                                      paste0(bundleName, bundle.annotation.suffix, '.json'))
            
            json = jsonlite::toJSON(bundleData, auto_unbox = TRUE, force = TRUE, pretty = TRUE)
            
            # use try mainly for permission problems on file system
            res=tryCatch(writeLines(json, annotFilePath), error=function(e) e)
            if(inherits(res,'error')){
              err=res
            }else{
              # annotation saved
              # set modified flag
              modified<<-TRUE
              
            }
            
            #### DBI ###
            # remove
            remove_bundleDBI(emuDBhandle, sessionName = bundleSession, name = bundleName)
            remove_bundleAnnotDBI(emuDBhandle, sessionName = bundleSession, bundleName = bundleName)
            # store
            # calculate MD5 sum of bundle annotJSON
            newMD5annotJSON = tools::md5sum(annotFilePath)
            names(newMD5annotJSON) = NULL
            
            bundleAnnotDFs = annotJSONcharToBundleAnnotDFs(as.character(json))
            add_bundleDBI(emuDBhandle, sessionName = bundleSession, name = bundleName, bundleAnnotDFs$annotates, bundleAnnotDFs$sampleRate, newMD5annotJSON)
            store_bundleAnnotDFsDBI(emuDBhandle, bundleAnnotDFs, sessionName = bundleSession, bundleName = bundleName)
            
            # rebuild redundant links & calculate posistions
            build_allRedundantLinks(emuDBhandle, sessionName = bundleSession, bundleName = bundleName)
            calculate_postionsOfLinks(emuDBhandle)
            
          }
        }
        if(is.null(err)){
          responseBundle=list(status=list(type='SUCCESS'),callbackID=jr$callbackID,responseContent='status',contentType='text/json')
        }else{
          m=err[['message']]
          cat('Error: ',m,"\n")
          responseBundle=list(status=list(type='ERROR',message=m),callbackID=jr[['callbackID']],responseContent='status',contentType='text/json')
        }
        # response object to JSON 
        responseBundleJSON=jsonlite::toJSON(responseBundle,auto_unbox=TRUE,force=TRUE,pretty=TRUE)
        # send response
        result=ws$send(responseBundleJSON)
        # restore warn level
        options(warn=warnOptionSave)
        
        # reset error
        err=NULL
        
      }else if(jr[['type']]=='DISCONNECTWARNING'){
        response=list(status=list(type='SUCCESS'),callbackID=jr[['callbackID']],responseContent='status',contentType='text/json')
        responseJSON=jsonlite::toJSON(response,auto_unbox=TRUE,force=TRUE,pretty=TRUE)
        result=ws$send(responseJSON)
        emuRserverRunning<<-FALSE
        ws$close()
        cat("emuR websocket service closed by EMU-webApp\n")
      }
    }
    ws$onMessage(serverReceive)
    ws$onClose(serverClosed)
  }
  
  app=list(call=httpRequest,onHeaders=onHeaders,onWSOpen=serverEstablished)
  sh=tryCatch(httpuv::startServer(host=host,port=port,app=app),error=function(e) e)
  if(inherits(sh,'error')){
    if(!is.null(getServerHandle())){
      cat("Trying to stop orphaned server (handle: ",getServerHandle(),")\n")
      httpuv::stopServer(getServerHandle())
      sh=tryCatch(httpuv::startServer(host=host,port=port,app=app),error=function(e) e)
      if(inherits(sh,'error')){
        stop("Error starting server (second try): ",sh,"\n")
      }
    }else{
      stop("Error starting server: ",sh,"\n")
    }
  }
  # store handle global for recovery after crash otr terminated R session
  setServerHandle(sh)
  cat("Navigate your browser to the EMU-webApp URL: http://ips-lmu.github.io/EMU-webApp/\n")
  cat("Server connection URL: ws://localhost:",port,"\n",sep='')
  cat("To stop the server press EMU-webApp 'clear' button or reload the page in your browser.\n")
  emuRserverRunning=TRUE
  if(length(autoOpenURL) != 0 && autoOpenURL != ""){
    # open browser with EMU-webApp
    utils::browseURL(autoOpenURL)
  }
  while(emuRserverRunning) {
    httpuv::service()
    Sys.sleep(0.01)
    
  }
  httpuv::stopServer(sh)
  # regular shutdown, remove handle 
  setServerHandle(NULL)
  if(debugLevel>0){
    cat("Closed emuR websocket HTTP service\n")
  }
  return(modified)
}

## searches for all tracks needed by the EMUwebApp and
## returns their ssffTrackDefinitions
get_ssffTracksUsedByDBconfig <- function(DBconfig){
  allTracks = NULL
  
  # anagestConfig ssffTracks
  for(ld in DBconfig$levelDefinitions){
    allTracks = c(allTracks, ld$anagestConfig$verticalPosSsffTrackName, ld$anagestConfig$velocitySsffTrackName)
  }
  
  for(p in DBconfig$EMUwebAppConfig$perspectives){
    # tracks in signalCanvases$order
    for(sco in p$signalCanvases$order){
      allTracks = c(allTracks, sco)
    }
    # tracks in twoDimCanvases$order
    for(tdco in p$twoDimCanvases$order){
      allTracks = c(allTracks, tdco)
    }
    
    # tracks in signalCanvases$assign
    for(sca in p$signalCanvases$assign){
      allTracks = c(allTracks, sca$ssffTrackName)
    }
    # tracks in p$twoDimCanvases$twoDimDrawingDefinitions
    for(tddd in p$twoDimCanvases$twoDimDrawingDefinitions){
      # dots
      for(dot in tddd$dots){
        allTracks = c(allTracks, dot$xSsffTrack, dot$ySsffTrack)
      }
    }
  }
  # uniq tracks
  allTracks = unique(allTracks)
  # remove OSCI and SPEC tracks
  allTracks = allTracks[allTracks != 'OSCI' & allTracks != 'SPEC']
  
  # get corresponding ssffTrackDefinitions
  allTrackDefs = list()
  for(std in DBconfig$ssffTrackDefinitions){
    if(std$name %in% allTracks){
      allTrackDefs[[length(allTrackDefs) + 1]] = std
    }
  }
  
  return(allTrackDefs)
}
