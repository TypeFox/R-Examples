


parseArgString <- function(string, delimiter=',', min.param.ct=2, max.param.ct=2, param.range=NULL){
## generic function for parsing delimited lists from BATCH mode argument strings.

  # count the delimiters
  re <- regexpr(delimiter, string)[[1]]
  delim.ct <- length(re)
  if(delim.ct >= max.param.ct)
    stop(paste('parameter count', delim.ct+1,'is too high for argument', string))
  if(re == -1 & min.param.ct != 1)
    stop(paste('you need a',delimiter,'delimited dimentions for this argument', string))

  p.vect <- strsplit(string, delimiter)[[1]]

  # check the range
  if(!is.null(param.range)){
    if(any(class(param.range[1]) == 'POSIXct'))    #if(class(param.range) == 'Date')
      p.vect <- as.POSIXct(p.vect)  #as.Date
    else
      p.vect <- as(p.vect, class(param.range))
    if(class(param.range) == 'factor')
      param.range <- as.character(param.range)
    if(class(param.range) == 'character'){
      if(!all(p.vect %in% param.range))
        stop(paste('not all of the elements of your parameter list', string, 'could be found in the allowed possiblities list', paste(param.range,collapse=',')))
    }else{
      if(!(param.range[1] <= range(p.vect)[1] & range(p.vect)[2] <= param.range[2]))
        stop(paste('the range you passed',string,'exceeds the allowed limits',paste(param.range,collapse=',')))
    }
  }
  return(p.vect)
}


.createBatchCommand <- function(cmdopts, script, logfile=paste(script,'.out',sep='')){
  cmd.rcmdbatch <- "R CMD BATCH --no-save --no-restore "
  cmd <- paste(cmd.rcmdbatch, cmdopts, script , logfile)
  return(cmd)
}


