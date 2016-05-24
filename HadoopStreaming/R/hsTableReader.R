`hsTableReader` <-
  function(file="",cols='character',chunkSize=-1,FUN=print, ignoreKey=TRUE,singleKey=TRUE, skip=0, sep='\t', keyCol='key', PFUN=NULL,carryMemLimit=512e6,carryMaxRows=Inf,stringsAsFactors=FALSE,debug=FALSE) {
    ## flush=TRUE  (allows comment, but precludes more than one record per line, which is good)
    if (skip > 0) {
      junk = scan(file,what='character',quiet=!debug,sep=sep,nlines=skip)
    }
    
    if (ignoreKey) {
      repeat {
        a = scan(file,what=cols,quiet=!debug,sep=sep,strip.white=TRUE,nlines=chunkSize,flush=TRUE)
        if ( length(a[[1]]) ==0 ) break
        FUN(data.frame(a,stringsAsFactors=stringsAsFactors))
      }
      return(invisible())
    }

    if (is.null(PFUN)) {
      ## Carryover Frame
      aCarry = data.frame()
      fileEmpty = TRUE
      repeat {
        if (object.size(aCarry)>carryMemLimit || nrow(aCarry) > carryMaxRows) {
          cat('In hsTableReader, aCarry has exceeded defined limits on memory and/or rows\n',file=stderr())
          cat('Key=',aCarry[1,keyCol],' MemoryUsed=',object.size(aCarry)/(2^20),'MB; NumRows=',nrow(aCarry),'\n',sep='',file=stderr())
          cat('Consider using higher values for carryMemLimit and carryMaxRows,\nOR use PFUN to handle incomplete keys.\n',file=stderr())
          ## Throw out the carryover data because getting too big
          aCarry=data.frame()
        }
        a = scan(file,what=cols,quiet=!debug,sep=sep,strip.white=TRUE,nlines=chunkSize,flush=TRUE)
        ## Memory Report
        if (debug) {
          cat('In hsTableReader, we have just scanned ',object.size(a)/(2^20),'MB. Current carry size is ',object.size(aCarry)/(2^20),'\n',file=stderr())
        }
        ## Done processing, because scan returned nothing
        if ( length(a[[keyCol]]) == 0 ) break
        fileEmpty = FALSE
        ## Prepend last carry to new data and convert scanned stuff to data.frame
        a = rbind(aCarry,data.frame(a,stringsAsFactors=stringsAsFactors))
        ##  Stick last user into aCarry
        lastKey = a[nrow(a),keyCol]
        firstRowOfLastKey = match(lastKey,a[,keyCol])
        aCarry = a[ firstRowOfLastKey:nrow(a) , ]
        if (!singleKey) {
          ## Process all complete keys at once
          FUN(a[1:(firstRowOfLastKey-1) , ])
          next
        }
        if (firstRowOfLastKey >= 2) {
          ## Process all complete keys, one at a time
          by(a[1:(firstRowOfLastKey-1),], a[1:(firstRowOfLastKey-1),keyCol], FUN)        
        }
      }
      if (!ignoreKey && !fileEmpty && nrow(aCarry)==0) stop ("empty aCarry at end -- this should never happen!!!")
      if (nrow(aCarry)>0) {
        FUN(aCarry)
      }
      return(invisible())
    }

    ## !is.null(PFUN)   here we handle partial keys in a streaming faction, rather than waiting for full key
    prevKey = NA
    repeat {
      a = scan(file,what=cols,quiet=!debug,sep=sep,strip.white=TRUE,nlines=chunkSize,flush=TRUE)
      if ( length(a[[keyCol]]) == 0 ) break
      a = data.frame(a,stringsAsFactors=stringsAsFactors)
      r = rle(a[,keyCol])
      numDistinctKeys = length(r$values)
      startKey = 1
      if (!is.na(prevKey)) {
        firstKey = r$values[1]
        if (prevKey == firstKey) {
          PFUN(a[1:r$lengths[1], ])
          if (numDistinctKeys == 1) {
            next
          }
          startKey = startKey + 1
        }
        ## We're guaranteed to be done with prevKey
        PFUN(data.frame())   #signals to PFUN that we're done the partial key
        prevKey = NA
      }
      firstRowOfLastKey = nrow(a) - r$lengths[numDistinctKeys] + 1
      startPos = if (startKey == 2) (1+r$lengths[1]) else 1
      if (startPos < firstRowOfLastKey) {
        if (!singleKey) {
          ## Process all complete keys at once
          FUN(a[startPos:(firstRowOfLastKey-1) , ])
        } else {
          ## Process all complete keys, one at a time
          for (keyNum in startKey:(numDistinctKeys-1)) {
            endPos = startPos+r$lengths[keyNum]-1
            FUN(a[startPos:endPos,])
            startPos = endPos+1
          }
          if (startPos != firstRowOfLastKey) stop("startPos != firstRowOfLastKey")
        }
      }
      prevKey =  a[ firstRowOfLastKey, keyCol ]
      PFUN( a[ firstRowOfLastKey:nrow(a) , ] )
    }
    if (!is.na(prevKey)) {
      PFUN(data.frame())  # signal to PFUN that we're done the last carry key
    }
    return(invisible())
  }
