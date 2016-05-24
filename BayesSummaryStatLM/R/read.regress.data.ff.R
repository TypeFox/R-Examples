read.regress.data.ff <- function(filename=NULL,predictor.cols=NA,response.col=NA,update.summaries=NULL
                              , fileEncoding = ""
                              , nrows = -1, first.rows = 1e5, next.rows = 1e5
                              , levels = NULL
                              , appendLevels = TRUE
                              , FUN = "read.table"
                              , transFUN = NULL
                              , asffdf_args = list()
                              , BATCHBYTES = getOption("ffbatchbytes")
                              , VERBOSE = FALSE
                              , header = FALSE
                              , sep = ","
                              , quote = "\"'"
                              , dec = "."
                              , numerals = c("allow.loss", "warn.loss", "no.loss")
                              , na.strings = "NA"
                              , colClasses = "numeric"
                              , skip = 0
                              , check.names = TRUE
                              , fill = TRUE
                              , strip.white = FALSE
                              , blank.lines.skip = TRUE
                              , comment.char = "#"
                              , allowEscapes = FALSE
                              , flush = FALSE
                              , skipNul = FALSE
)
{
  if(class(filename)=="character")
  {
    filelist <- list(filename)
  }
  else
  {
    if(class(filename)=="list")
      filelist <- filename
    else
      stop("read.regress.data: cannot interpret the filename argument")
  }
  regress.data <- NULL
  for (fname in filelist)
  {
    if((class(fname)=="character")&&(file.access(fname,4)==0))
    {
      finfo <- file.info(fname);
      if((!finfo$isdir)&&((finfo$mode&"444")>0)) # Check if the file is not a directory and has read permissions
      {
        data.ffdf<-read.table.ffdf(file=fname,
                              , fileEncoding = fileEncoding
                              , nrows = nrows, first.rows = first.rows, next.rows = next.rows
                              , levels = levels
                              , appendLevels = appendLevels
                              , FUN = FUN
                              , transFUN = transFUN
                              , asffdf_args = asffdf_args
                              , BATCHBYTES = BATCHBYTES
                              , VERBOSE = VERBOSE
                              , header = header
                              , sep = sep
                              , quote = quote
                              , dec = dec
                              , numerals = numerals
                              , na.strings = na.strings
                              , colClasses = colClasses
                              , skip = skip
                              , check.names = check.names
                              , fill = fill
                              , strip.white = strip.white
                              , blank.lines.skip = blank.lines.skip
                              , comment.char = comment.char
                              , allowEscapes = allowEscapes
                              , flush = flush
                              , skipNul = skipNul
        )
        
        numcols <- dim(data.ffdf)[2]
        N <- dim(data.ffdf)[1]
        
        if(((is.null(predictor.cols)||is.na(predictor.cols)))||!(is.numeric(predictor.cols))||!(length(predictor.cols)<numcols))
        {
          if((!is.na(predictor.cols))&&((!is.numeric(predictor.cols))||(length(predictor.cols)>numcols))) 
            warning('The parameter "predictor.cols" is either not properly formatted or does not match the data contained in file. The default settings are used.');
          predictor.cols <- 2:(numcols) # first column is assumed to be the response one
          response.col <- 1 #numcols
          numcols <- numcols-1
        }
        else
        {
          if(((is.null(response.col)||is.na(response.col)))||!(is.numeric(response.col))||!(length(response.col)=1))
          {
            response.col <- setdiff(1:numcols,predictor.cols)[1]
          }
          numcols <- length(predictor.cols)
        }
        # initialize the matrices for the loop below
        xtx = matrix(0.0,nrow = numcols, ncol = numcols)
        xty = matrix(0.0,nrow = numcols, ncol = 1)
        yty = 0.0;
        sum.x = rep(0.0,numcols);
        sum.y = 0.0;
        # go through the chunks of the data file
        for (chnk in chunk.ffdf(data.ffdf))
        {
          data.temp = data.ffdf[chnk,]
          x <- data.matrix(data.temp[predictor.cols])
          y <- data.matrix(data.temp[response.col])
          xtx <- xtx + t(x)%*%x 
          xty <- xty + t(x)%*%y 
          yty <- yty + t(y)%*%y
          sum.x <- sum.x + colSums(x)
          sum.y <- sum.y + sum(y)
        }
        
          numcols <- numcols+1
          xtx.ext              <- matrix(0,numcols,numcols)    
          xtx.ext[2:(numcols),2:(numcols)] <- xtx
          xtx.ext[1,1]         <- N
          xtx.ext[1,2:(numcols)]     <- t(sum.x)
          xtx.ext[2:(numcols),1]     <- sum.x
          
          xty.ext <- rep(0,(numcols))
          xty.ext[1]     <- sum.y
          xty.ext[2:(numcols)] <- xty 
          
          xtx <- xtx.ext
          xty <- xty.ext

        if(is.null(regress.data))
        {
            regress.data = list(xtx = xtx, xty = xty, yty = yty, numsamp.data = N)
        }
        else
          if(numcols == dim(regress.data$xtx)[1])
          {
            regress.data$xtx=regress.data$xtx+xtx
            regress.data$xty=regress.data$xty+xty
            regress.data$yty=regress.data$yty+yty
            regress.data$numsamp.data=regress.data$numsamp.data+N
          }
        else
          stop(paste("read.regress.data: The file \"",fname,"\" has data inconsistent with other files"));
      }
    else
    {
      warning(paste("read.regress.data: The file \"",fname,"\" cannot be read"));
    }
    }
    else
    {
      warning(paste("read.regress.data: The file \"",fname,"\" cannot be accessed"));
    }
  }
  if(!is.null(regress.data)&&(!(is.null(update.summaries)||is.na(update.summaries))))
  {
    if(all(!is.null(update.summaries[['xtx']]),is.matrix(update.summaries[['xtx']]),dim(update.summaries[['xtx']])==c(numcols,numcols),
           !is.null(update.summaries[['xty']]),is.numeric(update.summaries[['xty']]),length(update.summaries[['xty']])==numcols,
           !is.null(update.summaries[['yty']]),is.numeric(update.summaries[['yty']]),length(update.summaries[['yty']])==1,
           !is.null(update.summaries[['numsamp.data']]),is.numeric(update.summaries[['numsamp.data']]),length(update.summaries[['numsamp.data']])==1
    ))
    {
      regress.data$xtx=regress.data$xtx+update.summaries[['xtx']]
      regress.data$xty=regress.data$xty+update.summaries[['xty']]
      regress.data$yty=regress.data$yty+update.summaries[['yty']]
      regress.data$numsamp.data=regress.data$numsamp.data+update.summaries[['numsamp.data']]
    }
    else
      warning("read.regress.data.ff: incompatible update.summaries parameter is ignored.")
  }
  return(regress.data)
}