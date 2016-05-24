#! /usr/bin/env Rscript

## cat anna.R | ./hsWordCnt.R -m   | sort | ./hsWordCnt.R -r
library(HadoopStreaming)

## Additional command line arguments for this script (rest are default in hsCmdLineArgs)
spec = c('printDone','D',0,"logical","A flag to write DONE at the end.",FALSE)

opts = hsCmdLineArgs(spec, openConnections=T)

if (!opts$set) {
  quit(status=0)
}

mapperOutCols = c('word','cnt')
reducerOutCols = c('word','cnt')

if (opts$mapcols) {
  cat( paste(mapperOutCols,collapse=opts$outsep),'\n', file=opts$outcon )
} 

if (opts$reducecols) {
  cat( paste(reducerOutCols,collapse=opts$outsep),'\n', file=opts$outcon )
}

if (opts$mapper) {
  mapper <- function(d) {
    words <- strsplit(paste(d,collapse=' '),'[[:punct:][:space:]]+')[[1]] # split on punctuation and spaces
    words <- words[!(words=='')]  # get rid of empty words caused by whitespace at beginning of lines
    df = data.frame(word=words)
    df[,'cnt']=1
    hsWriteTable(df[,mapperOutCols],file=opts$outcon,sep=opts$outsep)
  }

  hsLineReader(opts$incon,chunkSize=opts$chunksize,FUN=mapper)

} else if (opts$reducer) {

  reducer <- function(d) {
    cat(d[1,'word'],sum(d$cnt),'\n',sep=opts$outsep)
  }
  cols=list(word='',cnt=0)  # define the column names and types (''-->string 0-->numeric)
  hsTableReader(opts$incon,cols,chunkSize=opts$chunksize,skip=opts$skip,sep=opts$insep,keyCol='word',singleKey=T, ignoreKey= F, FUN=reducer)
  if (opts$printDone) {
    cat("DONE\n");
  }
}

if (!is.na(opts$infile)) {
  close(opts$incon)
}

if (!is.na(opts$outfile)) {
  close(opts$outcon)
}



