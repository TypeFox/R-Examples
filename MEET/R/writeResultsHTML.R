writeResultsHTML <- function(resultsMEET, fileName='index.html') {
  require("MEET")
  require("seqinr")
  require("seqLogo")

  results <- resultsMEET$Results
  fileName <- 'index.html'
  iicc<-new.env()
  switch(resultsMEET$Summary$organism,
         "Rattus norvegicus"=	switch(resultsMEET$Summary$method,          	
                                     "Entropy"=data(RattusEntropy, envir=iicc),                                     
                                     "Divergence"=data(RattusDivergence, envir=iicc),                                     
                                     "Qresiduals"=data(RattusQresiduals, envir=iicc),
                                     stop("Method not included")),
         "Mus musculus"= 	switch(resultsMEET$Summary$method,                                 
                                 "Entropy"=data(MusEntropy, envir=iicc),                                     
                                 "Divergence"=data(MusDivergence, envir=iicc),                                     
                                 "Qresiduals"=data(MusQresiduals, envir=iicc),
                                 stop("Method not included")),
         "Drosophila melanogaster"= switch(resultsMEET$Summary$method,          	
                                           "Entropy"=data(DrosophilaEntropy, envir=iicc),                                     
                                           "Divergence"=data(DrosophilaDivergence, envir=iicc),                                     
                                           "Qresiduals"=data(DrosophilaQresiduals, envir=iicc),
                                           stop("Method not included")),
         "Homo sapiens"= 	switch(resultsMEET$Summary$method,          	
                                 "Entropy"=data(HomoEntropy, envir=iicc),                                     
                                 "Divergence"=data(HomoDivergence, envir=iicc),                                     
                                 "Qresiduals"=data(HomoQresiduals, envir=iicc),
                                 stop("Method not included")),
         stop("Organism not included"))
  
  filename <- "NoLogo"
  if (summary(paste("iicc", resultsMEET$Summary$nameTF,sep="")==ls(envir=iicc))[3]=="1") {    
    assign("TF",get(paste("iicc",resultsMEET$Summary$nameTF,sep=""),envir=iicc)$Transcriptionfactor)
    filename <- paste(fileName,".png",sep="")
    png(filename, bg = "transparent")
    seqLogo(con(TF,method="profile")/nrow(TF))
    dev.off()    
  } 

  pathMEET<-system.file(package="MEET")
  file.copy(paste(paste(pathMEET,"exec",sep="/"), "LogoMEET.png",sep="/"),"LogoMEET.png",overwrite=TRUE)
  
  sink(fileName)
  cat(c('<!doctype html><html><head><meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" /><title>WebMEET: A web of Multiple Detection of Transcription Factor Binding Sites</title>'), sep="\n")
  # webMEET.css
  cat(c('<style type="text/css">', readLines(paste(paste(pathMEET,"exec",sep="/"), "webmeet.css",sep="/"),encoding="UTF-8"), '</style>'), sep="\n")
  # jquery-ui-1.8.20.custom.css            
  cat(c('<style type="text/css">', readLines(paste(paste(pathMEET,"exec",sep="/"), "jquery-ui-1.8.20.custom.css",sep="/"),encoding="UTF-8"), '</style>'), sep="\n")
  # jquery.cluetip.css            paste(paste(pathMEET,"exec",sep="/"), "jquery.cluetip.css",sep="/")
  cat(c('<style type="text/css">', readLines(paste(paste(pathMEET,"exec",sep="/"), "jquery.cluetip.css",sep="/"),encoding="UTF-8"), '</style>'), sep="\n")
  # javascript
  jsArray <- 'var result = new Array('
  lines <- paste(paste(sapply(1:nrow(results), function(x) { 
      paste('{"Position":"', results[x,1], '","Value":"', results[x,2], '","Direction":"', results[x,3], '","Sequence":"', results[x,4], '","Organism":"', resultsMEET$Summary$organism, '","TF":"', resultsMEET$Summary$nameTF, '"}', sep='')
  })), collapse=",")
  jsArray <- paste(jsArray, lines, ');', sep='')
  cat(c('<script type="text/javascript">', jsArray, '</script>'), sep='')
  # jQuery.js             
  cat(c('<script type="text/javascript">', readLines(paste(paste(pathMEET,"exec",sep="/"), "jquery.js",sep="/"),encoding="UTF-8"), '</script>'), sep='\n')
  cat(c('<script type="text/javascript">', readLines(paste(paste(pathMEET,"exec",sep="/"), "jquery-ui-1.8.20.custom.min.js",sep="/"),encoding="UTF-8"), '</script>'), sep='\n')
  # jquery.cluetip.js                       
  cat(c('<script type="text/javascript">', readLines(paste(paste(pathMEET,"exec",sep="/"), "jquery.cluetip.js",sep="/"),encoding="UTF-8"), '</script>'), sep='\n')
  # webMEET.js                              
  cat(c('<script type="text/javascript">', readLines(paste(paste(pathMEET,"exec",sep="/"), "webmeet.js",sep="/"),encoding="UTF-8"),'</script>'), sep='')
  cat(c('</head><body><div id="main" class="inside"><div class="top_menu"><div class="logo"></div></div><div class="container">'), sep="\n")
  cat(c('<div class="tableContainer"><table id="result" class="tablesorter" cellspacing="1" cellpadding="0" border="0"><thead><tr><th class="header">Sequence</th><th class="header">p-value</th><th class="header">Position</th><th class="header">Direction</th></tr></thead><tbody class="scrollContent">'), sep="\n")
  cat(c('</tbody></table></div>'), sep="\n")
  cat(c('<span id="sequence_result" class="letter" style="display: inline;">', paste(sapply(1:length(resultsMEET$Summary$DNA[[1]]), function(i) { 
    letter <- resultsMEET$Summary$DNA[[1]][i]
    if (i %in% results[,1]) {
      pos <- match(i, results[,1])
      letter <- paste("<span id=\"", results[pos,4], results[pos,1], results[pos,3], "\" class=\"mark\" title=\"|<span class=\'tit\'>Organism:</span> <span class=\'res\'>", resultsMEET$Summary$organism, "</span>|<span class=\'tit\'>TF:</span> <span class=\'res\'>", resultsMEET$Summary$nameTF, "</span>|<span class=\'tit\'>Sequence:</span> <span class=\'res\'>", results[pos,4], "</span>|<span class=\'tit\'>p-value:</span> <span class=\'res\'>", results[pos, 2], "</span>|<span class=\'tit\'>Position:</span> <span class=\'res\'>", results[pos, 1], "</span>|<span class=\'tit\'>Direction:</span> <span class=\'res\'>", results[pos, 3], "\"</span>", letter, "</span>", sep='')
    }
    if (i %% 60 == 0) {
      paste(letter, '\n', sep='')
    } else {
      letter
    }
  }), collapse=''), '</span>'), sep='\n')
  cat(c('</div></div><div id="dialog" title="Basic dialog"><img id="consensus" src="', filename, '" style="width:350px" /></div></body></html>'))
  sink()
 # browseURL(paste('file:///', getwd(), fileName, sep='/'))
}