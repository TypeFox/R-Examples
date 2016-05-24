#' nat.utils: Defunct functions
#' 
#' @description These functions have been retired from \code{nat.utils}
#' @param from Source file
#' @param to (New) target hardlink file to create
#' @export
#' @name nat.utils-defunct
#' @aliases file.hardlink
#' @seealso \code{\link{file.link}}
file.hardlink=function(from, to){
  .Defunct('file.link', package = 'nat.utils')
}

#' Swap names of two files (by renaming first to a temporary file)
#' @param f1,f2 Paths to files
#' @return logical indicating success
#' @author jefferis
#' @export
#' @seealso \code{\link{file.rename}}
file.swap<-function(f1,f2){
  # quick function to swap filenames 
  if(length(f1)>1 || length(f2)>1) return(mapply(file.swap,f1,f2))
  
  if(!all(file.exists(f1),file.exists(f2))) stop("f1 and f2 must exist")
  
  tmpfile=tempfile(basename(f1),dirname(f1))
  rval=file.rename(from=f1,to=tmpfile) && file.rename(from=f2,to=f1) && file.rename(from=tmpfile,to=f2)
  return(rval)
}

#' Remove common part of two paths, leaving relative path
#' @param path Paths to make relative
#' @param stempath Root to which \code{path} will be made relative
#' @param StopIfNoCommonPath Error if no path in common
#' @return Character vector containing relative path
#' @author jefferis
#' @export
#' @seealso \code{\link{path.expand}}, \code{\link{normalizePath}}
#' @family path_utils
#' @examples
#' path = "/Volumes/JData/JPeople/Sebastian/images"
#' abs2rel(path,'/Volumes/JData')
abs2rel<-function(path,stempath=getwd(),StopIfNoCommonPath=FALSE){
  if(length(stempath)>1)
    stop("only 1 stempath allowed!")
  path=path.expand(path)
  stempath=path.expand(stempath)
  ncsp=nchar(stempath)
  fsep=.Platform$file.sep
  if(substr(stempath,ncsp,ncsp)!=fsep){
    stempath=paste(stempath,fsep,sep='')
  }
  
  relpath=sub(stempath,"",path,fixed=TRUE)
  
  warnorstopfun=if(StopIfNoCommonPath) stop else warning
  
  mismatches=which(relpath==path)
  if(length(mismatches)){
    warnorstopfun("stempath: ",stempath," is not present in path(s): ",
                  paste(path[mismatches], collapse=":"))
  }
  relpath
}

#' Find common prefix of two or more (normalised) file paths
#' 
#' @param paths Character vector of file paths
#' @param normalise Whether to normalise \code{paths} (with 
#'   \code{\link{normalizePath}}, default \code{FALSE})
#' @param fsep Optional path separator (defaults to \code{.Platform$file.sep})
#' @return Character vector of common prefix, \code{""} when there is no common 
#'   prefix, or the original value of \code{paths} when fewer than 2 paths were 
#'   supplied.
#' @export
#' @details Note that for absolute paths, the common prefix will be returned 
#'   e.g. \code{common_path(c("/a","/b"))} is \code{"/"}
#'   
#'   Note that \code{\link{normalizePath}} 1) operates according to the 
#'   conventions of the current runtime platform 2) is called with 
#'   \code{winslash=.Platform$file.sep} which means that normalised paths will 
#'   eventually end up separated by "\" by default on Windows rather than by 
#'   "//", which is \code{normalizePath}'s standard behaviour.
#' @seealso \code{\link{normalizePath}}
#' @family path_utils
#' @examples
#' common_path(c("/a","/b"))
#' common_path(c("/a/b/","/a/b"))
#' common_path(c("/a/b/d","/a/b/c/d"))
#' common_path(c("/a/b/d","/b/c/d"))
#' common_path(c("a","b"))
#' common_path(c("","/a"))
#' common_path(c("~","~/"))
#' common_path(c("~/a/b/d","~/a/b/c/d"), normalise = FALSE)
#' common_path(c("~","~/"), normalise = FALSE)
common_path<-function(paths, normalise=FALSE, fsep=.Platform$file.sep) {
  if(normalise)
    paths=normalizePath(paths, mustWork = F, winslash = fsep)
  
  if(length(paths)<2) 
    return(paths)

  path_chunks=lapply(paths, split_path, fsep=fsep, include.fseps=TRUE, 
                     omit.duplicate.fseps=TRUE)
  maxlen=max(sapply(path_chunks, length))
  # pad cols with NAs to same length
  path_chunks=lapply(path_chunks, function(x) {length(x)=maxlen;x})
  # make a char matrix with 1 col per path and one row per chunk
  m=do.call(cbind, path_chunks)
  # figure out who is equal row by row
  num_uniq_values=apply(m, 1, function(x) length(unique(x)))
  # establish which chunks are different
  diff_chunks=which(num_uniq_values!=1)
  first_diff_chunk=min(c(diff_chunks, maxlen+1))
  # paste the chunks upt to the first different fragment
  paste(m[seq_len(first_diff_chunk-1), 1], collapse="")
}


#' Split file path into individual components (optionally including separators)
#' 
#' @param path A path with directories separated by \code{fsep}s.
#' @param include.fseps Whether to include the separators in the returned 
#'   character vector (default \code{FALSE})
#' @param omit.duplicate.fseps Whether to omit duplicate file separators if 
#'   \code{include.fseps=TRUE} (default \code{FALSE}).
#' @param fsep The path separator (default to \code{.Platform$file.sep})
#'   
#' @return A character vector with one element for each component in the path 
#'   (including path separators if \code{include.fseps=TRUE}).
#' @export
#' @family path_utils
#' @seealso \code{\link{file.path}}
#' @examples
#' split_path("/a/b/c")
#' split_path("a/b/c")
#' parts=split_path("/a/b/c", include.fseps=TRUE)
#' # join parts back up again
#' paste(parts, collapse = "")
#' split_path("a/b//c", include.fseps=TRUE, omit.duplicate.fseps=TRUE)
#' # Windows style
#' split_path("C:\\a\\b\\c", fsep="\\")
split_path<-function(path, include.fseps=FALSE, omit.duplicate.fseps=FALSE,
                     fsep=.Platform$file.sep) {
  if(length(path)>1) 
    stop("I can only handle one path")
  if(nchar(fsep)>1)
    stop("fsep must consist of one character only!")
  
  # nb c() clears attributes
  seps=c(gregexpr(fsep, path, fixed = T)[[1]])
  # no match
  if(seps[1]<0) return(path)
  # "/a/b/c" -> "/" "a" "/" "b" "c"
  # "job/cat/" -> "job" "/" "cat" "/"
  
  # add fake separator in last position to simplify loop below
  seps=c(seps, nchar(path)+1)
  chunks=character()
  p=1
  while(p<=nchar(path)){
    if(substr(path,p,p)==fsep) {
      if(include.fseps) {
        if(omit.duplicate.fseps && isTRUE(chunks[length(chunks)]==fsep)) {
          # omit this chunk
        } else chunks=c(chunks, fsep)
      }
      p=p+1
    } else {
      # grab everything up to next separator
      nextsep=seps[seps>p][1]
      chunk=substr(path, p, nextsep-1)
      p=p+nchar(chunk)
      chunks=c(chunks,chunk)
    }
  }
  chunks
}

#' Use unix touch utility to change file's timestamp
#' 
#' If neither a time or a reference file is provided then the current time is 
#' used. If the file does not already exist, it is created unless Create=FALSE.  
#' @param file Path to file to modify
#' @param time Absolute time in POSIXct format
#' @param reference Path to a reference file
#' @param timestoupdate "access" or "modification" (default both)
#' @param Create Logical indicating whether to create file (default TRUE)
#' @return TRUE or FALSE according to success
#' @author jefferis
#' @export
touch<-function(file,time,reference,timestoupdate=c("access","modification"),
    Create=TRUE){
  if(.Platform$OS.type!="unix")
    stop("touch relies on the existence of a system touch command")

  if(!Create && !file.exists(file)) stop("Create=F and ",file," does not exist") 
  if(!missing(time) && !missing(reference))
    stop("Please supply either a time or a reference file but not both")
  args=paste("-",substr(timestoupdate,1,1),sep="")
  if(!missing(time)){
    if(!is.character(time)) time=strftime(time,"%Y%m%d%H%M.%S")
    args=c(args,"-t",time)
  } else if(missing(reference)) {
    # use current time
  } else {
    # use reference file to supply time
    if(!file.exists(reference)) stop("reference file: ",reference," missing")
    args=c(args,"-r",shQuote(reference))
  }
  
  cmd=paste("touch",paste(args,collapse=" "),shQuote(file))
  return(system(cmd)==0)
}
