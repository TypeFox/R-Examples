#' Return information about a zip archive using system unzip command
#' 
#' @details Uses system unzip command.
#' @param f Path to one (or more) files
#' @return dataframe of information
#' @family ziputils
#' @author jefferis
#' @export
#' @seealso \code{\link{zip}}
#' @importFrom utils read.table
zipinfo<-function(f){
  if(length(f)>1) return(sapply(f,zipinfo))
  results=system2(command=unzip(T), args=c("-lv", shQuote(f)), stdout=TRUE)
  
  dash_lines=grep('^-----',results)
  if(length(dash_lines)!=2) stop("Unable to parse zip information for file ",f)
  result_lines=seq(from=dash_lines[1]+1,to=dash_lines[2]-1)
  
  selresults=results[c(2,result_lines)]
  nameColPos=regexpr('Name',results[2])
  firstcols=substr(selresults,1,nameColPos-1)
  
  # remove initial spaces
  firstcols=sub("^[ ]+",'',firstcols)
  # convert intervening spaces to tabs
  firstcols=gsub("[ ]+",'\t',firstcols)
  df=read.table(text=firstcols,header=T,colClasses=c(Method='factor'),as.is=TRUE)
  
  # now handle file names
  max_width=max(sapply(results,nchar))
  df$Name=substr(selresults[-1],nameColPos,max_width)
  df
}

#' Verify integrity of one or more zip files
#' 
#' @details Uses system unzip command.
#' @param f Path to one (or more) files
#' @param Verbose Whether to be Verbose (default FALSE)
#' @return TRUE when file OK, FALSE otherwise
#' @author jefferis
#' @export
#' @family ziputils
zipok<-function(f,Verbose=FALSE){
  if(length(f)>1) return(sapply(f, zipinfo, Verbose=Verbose))
  result=system2(command=unzip(T), args=c("-t", shQuote(f)), 
                 stdout=ifelse(Verbose, "", FALSE))
  result==0
}

unzip<-function(mustWork=FALSE) {
  if(is.null(w<-getOption("nat.utils.unzip"))){
    w=Sys.which("unzip")[[1]]
    options(nat.utils.unzip=w)
  }
  if(mustWork && !nzchar(w))
    stop("Cannot find system unzip command!")
  return(w)
}
