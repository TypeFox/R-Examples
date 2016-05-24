## Simply make this accessible here as a convenience to the user

#' @export
#' @importFrom utils download.file

lookup.xport <- function(file)
  {
    fname <- file

    if(length(grep('http://', file))>0 || length(grep('ftp://', file))>0 )
      {
        scat("Downloading file...")        
        tf <- tempfile()
        download.file(file, tf, mode='wb', quiet=TRUE)
        file <- tf
      }

    ret <- lookup.xport.inner(file)
    attr(ret, "call") <- match.call()
    attr(ret, "file") <- fname
    class(ret) <- c("lookup.xport", "list")
    ret
  }

#' @export

print.lookup.xport <- function(x, ...)
  {
    Sinfo <- summary(x, ...)
    print(Sinfo)
  }

#' @export

summary.lookup.xport <- function(object, ...)
  {
    subFun <- function(XX)
      {
        df <- object[[XX]]
        ret <- as.data.frame(df[c(
                                  "name", "type",
                                  "format", "flength", "fdigits",
                                  "iformat", "iflength", "ifdigits",
                                  "label"
                                  )
                                ])
        if(nrow(ret)==0) ret[1,] <- NA # ensure at least one row
        cbind(dataset=XX, ret, nobs=df$length)
      }
    
    dFrames <- lapply( names(object), subFun )
    singleFrame <- do.call("rbind", dFrames)
    rownames(singleFrame) <- paste(singleFrame$dataset, singleFrame$name, sep=".")

    attr(singleFrame, "call") <- attr(object, "call")
    attr(singleFrame, "file") <- attr(object, "file")    
    class(singleFrame) <- c("summary.lookup.xport","data.frame")
    
    singleFrame
  }

#' @export

print.summary.lookup.xport <- function(x, ...)
{
  cat("\n")
  cat("SAS xport file\n")
  cat("--------------\n");
  cat("Filename: `", attr(x,"file"), "'\n", sep="")
  cat("\n")
  for(dSetName in levels(x$dataset))
    {
      cat("Variables in data set `", dSetName, "':\n", sep="")
      print(as.data.frame(x)[x$dataset==dSetName,], row.names=FALSE)
      cat("\n")
    }
}

