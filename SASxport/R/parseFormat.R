## Convert SAS format specification string into format name, length, and digits
parseFormat <- function(format)
  {
    retval <- list("name"="", "len"=0, "digits"=0)


    if( !is.null(format) && (length(format)==1) && (format > "") )
              {
                index <- regexpr("[0-9]+", format)
                if(index==-1)
                  {
                    retval$name   <- format
                    retval$len    <- 0
                    retval$digits <- 0
                  }
                else
                  {
                    retval$name <- substr(format,0,index-1)[1]

                    lenStr <- substr(format, index, nchar(format))

                    index <- regexpr("\\.", lenStr)
                    if(index==-1)
                      {
                        retval$len    <- as.numeric(lenStr)
                        retval$digits <- 0
                      }
                    else
                      {
                        retval$len     <- as.numeric(substr(lenStr, 0, index-1))
                        retval$digits  <- as.numeric(substr(lenStr, index+1, nchar(lenStr)))
                      }
                  }

                if(is.na(retval$len)) retval$len <- 0
                if(is.na(retval$digits)) retval$digits <- 0

              }

    return(retval)
  }
