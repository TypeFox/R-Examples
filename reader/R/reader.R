###NAMESPACE ADDITIONS###
# Depends: R (>= 2.10), utils, NCmisc
# Imports: grDevices, graphics, stats
# Suggests:
# importFrom(stats, rnorm)
# importFrom(graphics, hist, plot)
# importFrom(grDevices, colors) 
# import(utils, NCmisc)
###END NAMESPACE###

#' Find which column in a dataframe contains a specified set of values.
#' 
#' Starting with a list of ids, each column is searched. The column with
#'  the highest non-zero percentage matching is assumed to correspond
#'  to the id list. The search terminates early if a perfect match is 
#'  found. Useful for assembling annotation from multiple sources.
#'
#' @param frame a data.frame, or similarly 2 dimensional object which 
#'  might contain ids
#' @param ids a vector of IDs/value that might be found in at least 
#'  1 column of frame
#' @param ret specify what should be returned, see values
#' @return ret can specify a list returning, 'col': the column number 
#'  (col=0 for rownames) with the best match; 'maxpc': the percentage
#'  of ids found in the best matching column; 'index': the matching vector
#'  that maps the frame rows onto ids; 'results': the (sub)set of ids
#'  found in frame. NAs given for ids not found
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' new.frame <- data.frame(day=c("M","T","W"),time=c(9,12,3),staff=c("Mary","Jane","John"))
#' staff.ids <- c("Mark","Jane","John","Andrew","Sally","Mary")
#' new.frame; staff.ids; find.id.col(new.frame,staff.ids)
find.id.col <- function(frame,ids,ret=c("col","maxpc","index","result"))
{
  # looks across each column of a dataframe 'frame' and selects the
  # column that has the most matches to the vector 'ids'
  # can return any/all of:
  #   the number of best match column, the match %, index numbers, matches
  opts <- c("col","maxpc","index","result")
  if(is.null(dim(frame))) { warning("object 'frame' has no dimensions, not suitable for find.id.col()"); return(NULL) }
  targ.cols <- ncol(frame)+1 # +1 if for looking at 'rownames'
  if(targ.cols==1) { warning("data.frame had no columns"); return(NULL) }
  pcs <- numeric(targ.cols)
  num.ids <- length(ids); if(num.ids <1) { warning("at least 1 ID must be entered") ; return(NULL) }
  ids <- paste(ids)
  if(!is.character(ids)) { stop("ids must be coercible to character type") }
  coln <- 1; best <- NA
  for (cc in 1:targ.cols)
  {
    if(cc==targ.cols)
    {
      # last 'cc' is actually to test the rownames of the frame
      if(!is.null(rownames(frame))) { 
        if(all(rownames(frame)==paste(1:nrow(frame)))) {
          posit <- rep(NA,length(ids)) # force NAs if rownames are just column number
          #  warning("rownames were just column numbers so NA match returned")
        } else {
          posit <- match(ids,rownames(frame)) 
        }
      } else { break }
    } else {
      posit <- match(ids,frame[,cc])
    }
    pcs[cc] <- length(which(!is.na(posit)))/num.ids
    if(pcs[cc]>max(pcs[-cc])) { best <- posit ; coln <- cc }
    if(pcs[cc]==1) { break } # exit if found a perfect match
  }
  maxpc <- max(pcs)
  if(coln==targ.cols) {
    result <- rownames(frame)[best]
    coln <- 0
  } else {
    result <- frame[best,coln]
  }
  out <- list(coln,maxpc,best,result)
  names(out) <- opts
  for (cc in length(opts):1)
  {
    if (!(opts[cc] %in% ret)) { out[[cc]] <- NULL }
  }
  return(out)
}


#' returns a dataframe if 'unknown.data' can in anyway relate to such:
#'
#' it can be:
#' - dataframe, matrix, big.matrix, sub.big.matrix, big.matrix.descriptor,
#' a bigmatrix description file, an RData file containing one of these 
#' objects, the name of a text or RData file, a named vector (names 
#' become rownames), or a list containing a matrix or dataframe. 
#' Using this within functions allows flexibility in 
#' specification of a datasource
#'
#' @param unknown.data something that is or can refer to a 2d dataset
#' @param too.big max size in GB, to prevent unintended conversion to 
#'  matrix of a very large big.matrix object.
#' @return returns a data.frame regardless of the original object type
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{force.vec}}
#' @examples
#' # create a matrix, binary file, text file, big.matrix.descriptor
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' test.files <- c("temp.rda","temp.txt")
#' mymat <- matrix(rnorm(100),nrow=10)
#' # not run yet # require(bigmemory)
#' save(mymat,file=test.files[1])
#' write.table(mymat,file=test.files[2],col.names=FALSE,row.names=FALSE)
#' test.frames <- list(mymat = mymat,
#'  myrda = test.files[1], mytxt = test.files[2] )
#'  # not run yet #: ,mybig = describe(as.big.matrix(mymat)) )
#' sapply(sapply(test.frames,is),"[",1)
#' # run the function on each, reporting specs of the object returned
#' for (cc in 1:length(test.frames)) {
#'   the.frame <- force.frame(test.frames[[cc]])
#'   cat(names(test.frames)[cc],": dim() => ",
#'       paste(dim(the.frame),collapse=","),
#'       "; is() => ",is(the.frame)[1],"\n",sep="")
#' }
#' unlink(test.files)
#' setwd(orig.dir) # reset working dir to original
force.frame <- function(unknown.data,too.big=10^7)
{
  # returns a dataframe if 'unknown.data' can in anyway relate to such:
  # it can be:
  # - dataframe, matrix, big.matrix, sub.big.matrix, big.matrix.descriptor,
  # a bigmatrix description file, an RData file containing one of these objects,
  # the name of a text or RData file, a named vector (names become rownames),
  # or a list containing a matrix or dataframe. Using this in functions allows
  # flexibility in specification of a datasource
  uis <- is(unknown.data)
  if(uis[1]=="character" & length(unknown.data)==1)
  {
    # filename?
    if(file.exists(unknown.data)) { 
      out.data <- reader(unknown.data)
    } else {
      stop("Error: argument seemed to be a filename (char, length 1) but file did not exist")
    }
  } else {
    if(uis[1] %in% c("matrix","data.frame")) {
      out.data <- unknown.data
    } else {
      if(uis[1]=="list") {
        types <- sapply(unknown.data,is)
        wty <- which(types %in% c("matrix","data.frame"))[1]
        if(length(wty)==1) {
          out.data <- unknown.data[[wty]]
        } else {
          stop("Error: object was list and no elements were a matrix or data.frame")
        }
      } else {
        if(is.null(dim(unknown.data)) & !is.null(names(unknown.data)))
        {
          cat(" converting named vector into dataframe\n")
          out.data <- as.matrix(unknown.data)
        } else {
          if(length(grep("big.matrix",uis))>0) {
            if(!exists("get.big.matrix",mode="function")) { 
              warning("'bigpc' package missing can't read big.matrix"); return(NULL) }
            unknown.data <- do.call("get.big.matrix",list(bigMat=unknown.data)) # once bigpc  is submitted replace with actual fn
            if(estimate.memory(unknown.data) <= too.big) {
              cat(" converting big.matrix object into a dataframe\n")
              out.data <- as.matrix(unknown.data) 
            } else {
              stop(paste("Error: big matrix object was too big to convert to a dataframe [>",
                         too.big,"cells]\n"))
            }
          } else {
            warning(paste("trying to convert object type",uis[1],"into a dataframe"),
                    " - result may be unpredictable")
            out.data <- as.df(unknown.data)
            #print(head(out.data))
            return(out.data)
          }
        }
      }
    }
  }
  return(as.df(out.data))
}


#' Internal function to extract the best guess at an ID vector from
#' a dataframe, without knowing any ids (set most.unique=FALSE to select col 1)
#' default to select the most unique column from the first 100 cols
vec.extract.mat <- function(X,most.unique=TRUE,max.col=100) {
  num.unique <- function(X) { length(unique(X)) }
  if(!is.null(rownames(X))) {
    if(!all(rownames(X)==paste(1:nrow(X)))) {
      return(rownames(X))
    }
  } 
  if(most.unique) {
    unz <- sapply(X[1:min(max.col,ncol(X))],num.unique)
    index <- which(unz==max(unz,na.rm=TRUE))[1]
    if(length(index)==0) { index <- 1 }
  } else {
    index <- 1
  }
  return(X[,index])
}



#' returns a vector if 'unknown.data' can in anyway relate to such:
#'
#' if the name of a file with a vector or vector, then reads the file,
#' if a matrix or dataframe, then preferentially return rownames, 
#'  otherwise return first column - designed to search for IDs.
#' Using this within functions allows flexibility in 
#' the specification of a datasource for vectors
#'
#' @param unknown.data something that is or can refer to a 2d dataset
#' @param most.unique if TRUE, select most unique column if a unknown.data
#'  is a matrix, else select the first column
#' @param dir if unknown.data is a file name, specifies directory(s) to
#'  look for the file
#' @param warn whether to display a warning if unknown.data is a matrix
#' @return returns a vector regardless of the original object type
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{force.frame}}
#' @examples
#' # create a matrix, binary file, and simple vector
#' my.ids <- paste("ID",1:4,sep="")
#' my.dat <- sample(2,4,replace=TRUE)
#' test.files <- c("temp.rda")
#' mymat <- cbind(my.ids,my.dat)
#' save(mymat,file=test.files[1])
#' test.vecs <- list(myvec = my.ids,
#'  myrda = test.files[1],mymat=mymat)
#' # show dimensions of each test object
#' sapply(test.vecs,function(x) {  if(is.null(dim(x))){ length(x)} else {dim(x)}})
#' # run the function on each, reporting specs of the object returned
#' for (cc in 1:3) {
#'   the.vec <- force.vec(test.vecs[[cc]])
#'   cat(names(test.vecs)[cc],": length() => ",
#'       length(the.vec),"; is() => ",is(the.vec)[1],"\n",sep="")
#' }
#' unlink(test.files)
force.vec <- function(unknown.data,most.unique=TRUE,dir=NULL, warn=FALSE)
{
  # returns a vector if 'unknown.data' can in anyway relate to such:
  # if the name of a file with a vector or vector, then returns that
  # if a matrix or dataframe, then preferentially return rownames, otherwise
  # return first column
  # designed to search for IDs basically
  # 
  uis <- is(unknown.data)
  if(is.null(dim(unknown.data))) {
    out.data <- unknown.data
    if(!"vector" %in% uis) {
      warning("unknown dimensionless non-vector datatype!")
    } else {
      if(uis[1]=="character" & length(unknown.data)==1) {
        if(is.file(unknown.data,dir)) { 
          out.data <- reader(find.file(unknown.data,dir))
          #print(head(out.data))
          out.data <- force.vec(out.data)
        } else {
          warning("argument length 1 might be a filename (char, length 1) but file did not exist")
          cat("arg: \n"); print(unknown.data)
        }
      } 
    }
    if("list" %in% uis) { out.data <- unlist(out.data) }
  } else {
    out.data <- force.frame(unknown.data)
    if(warn) {
      warning("input seems to be dim>1, rather than a vector, using 'vec.extract.mat' to find best vector")
    }
    out.data <- vec.extract.mat(out.data,most.unique)
  }
  return(as.vector(out.data))
}



#' Change column name in different form to desired form.
#' 
#' Searches for possible equivalents for a desired column in a dataframe
#' and replaces first name match with desired name. Useful when parsing
#' different annotation files which may have standard columns with slightly
#' different names, e.g, Gender=SEX=sex=M/F, or ID=id=ids=samples=subjectID
#'
#' @param frame a dataframe or matrix with column names
#' @param desired the column name wanted
#' @param testfor possible alternate forms of the desired column name
#' @param ignore.case whether to ignore the upper/lower case of the column names
#' @return returns the original dataframe with the target column renamed
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' df <- data.frame(Sex=c("M","F","F"),time=c(9,12,3),ID=c("ID3121","ID3122","ID2124"))
#' # standard example
#' new.df <- column.salvage(df,"sex",c("gender","sex","M/F")); df; new.df
#' # exact column already present so no change
#' new.df <- column.salvage(df,"ID",c("ID","id","ids","samples","subjectID")); df; new.df
#' # ignore case==TRUE potentially results in not finding desired column:
#' new.df <- column.salvage(df,"sex",c("gender","sex","M/F"),ignore.case=FALSE); df; new.df
column.salvage <- function(frame,desired,testfor, ignore.case=TRUE) 
{
  ## attempt to find predictable misnaming of column (contents of 'testfor')
  # and change name of the first detected misname in the 'testfor' list to 'desired'
  # e.g, want desired column 'GRP' and look for misnamings: 'group', 'Grp', 'grp', etc
  if(is.null(colnames(frame))) { warning("frame had no column names"); return(frame) }
  if(!all(desired %in% colnames(frame))) {
    # ^ ie, if 'desired' already present, do nothing
    tf <- testfor; cf <- colnames(frame)
    if(ignore.case) { tf <- c(tolower(desired),tolower(tf)); cf <- tolower(cf) } 
    if(any(tf %in% cf)) {
      colnames(frame)[(narm(match(tf,cf)))[1]] <- desired
    } else {
      warning("couldn't find any columns: ",paste(testfor,collapse=", ")," to change to '",desired,"' in frame")
    }
  }
  return(frame)
}






#' Flexibly load from a text or binary file, accepts multiple file formats.
#' 
#' Uses file extension to distinguish between binary, csv or other 
#' text formats. Then tries to automatically determine other parameters
#' necessary to read the file. Will attempt to detect the delimiter,
#' and detect whether there is a heading/column names, and whether 
#' the first column should be rownames, or left as a data column.
#' Internal calls to standard file reading functions use 
#' 'stringsAsFactors=FALSE'.
#'
#' @param fn filename (with or without path if dir is specified)
#' @param dir optional directory if separate path/filename is preferred
#' @param want.type if loading a binary file with multiple objects, specify
#'  here the is() type of object you are trying to load
#' @param def the default delimiter to try first
#' @param force.read attempt to read the file even if the file type looks unsupported
#' @param header presence of a header should be autodetected, but can specify header status 
#' if you don't trust the autodetection
#' @param h.test.p p value to discriminate between number of characters in a column name versus
#'  a column value (sensitivity parameter for automatic header detection)
#' @param quiet run without messages and warnings
#' @param treatas a standard file extension, e.g, 'txt', to treat file as
#' @param override assume first col is rownames, regardless of heuristic
#' @param more.types optionally add more file types which are read as text
#' @param auto.vec if the file seems to only have a single column, automatically
#'  return the result as a vector rather than a dataframe with 1 column
#' @param one.byte logical parameter, passed to 'get.delim', whether to look for only 1-byte
#'  delimiters, to also search for 'whitespace' which is a multibyte (wildcard) delimiter type. 
#'  Use one.byte = FALSE, to read fixed width files, e.g, many plink files.
#' @param ... further arguments to the function used by 'reader' to parse the file,
#'  e.g, depending on file.type, can be read.table(), read.delim(), read.csv().
#' @return returns the most appropriate object depending on the file type,
#'  which is usually a data.frame except for binary files
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' # create some datasets
#' df <- data.frame(ID=paste("ID",101:110,sep=""),
#'   scores=sample(70,10,TRUE)+30,age=sample(7,10,TRUE)+11)
#' DNA <- apply(matrix(c("A","C","G","T")[sample(4,100,TRUE)],nrow=10),
#'                                                 1,paste,collapse="")
#' fix.wid <- c("    MyVal    Results        Check",
#'   "    0.234      42344          yes",
#'   "    0.334        351          yes","    0.224         46           no",
#'   "    0.214     445391          yes")
#' # save data to various file formats
#' test.files <- c("temp.txt","temp2.txt","temp3.csv",
#'                               "temp4.rda","temp5.fasta","temp6.txt")
#' write.table(df,file=test.files[1],col.names=FALSE,row.names=FALSE,sep="|",quote=TRUE)
#' write.table(df,file=test.files[2],col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)
#' write.csv(df,file=test.files[3])
#' save(df,file=test.files[4])
#' writeLines(DNA,con=test.files[5])
#' writeLines(fix.wid,con=test.files[6])
#' # use the same reader() function call to read in each file
#' for(cc in 1:length(test.files)) {
#'   cat(test.files[cc],"\n")
#'   myobj <- reader(test.files[cc])  # add 'quiet=FALSE' to see some working
#'   print(myobj); cat("\n\n")
#' }
#' # inspect files before deleting if desired
#' unlink(test.files) 
#' # myobj <- reader(file.choose()); myobj # run this to attempt opening a file
#' setwd(orig.dir) # reset working directory to original
reader <- function(fn,dir="",want.type=NULL,def="\t",force.read=TRUE,header=NA,h.test.p=0.05,
                   quiet=TRUE,treatas=NULL,override=FALSE,more.types=NULL,
                   auto.vec=TRUE,one.byte=TRUE,...)
{
  # try to read in data from many types of datafile
  typ <- classify.ext(fn,more.txt=more.types)
  #print(typ)
  types <- c("BIN","CSV","TXT","OTH")
  treatas <- classify.ext(treatas,more.txt=more.types)
  # get file extension
  fsuf <- get.ext(fn)
  # construct the full path
  full.path <- cat.path(dir,fn)
  if(!file.exists(full.path)) { stop("file 'fn' did not exist") }
  # select file type to open as
  # typ <- types[types %in% fsuf][1]
  if(!is.null(treatas) )
  {
    if(treatas[1] %in% types) {
      cat(paste(" will assume file is",treatas,"\n"))
      typ <- treatas[1]
    } 
  }
  if(typ=="OTH") {
    if(force.read) {
      if(!quiet) { warning("looks like a non-supported file type") }
      typ <- types[3]
    } else {
      warning("looks like a non-supported file type.",
              "\nadd suffix to 'more.types' if it's a text file",
              "\nchange file suffix to .rda or .csv if it's binary or csv")
      cat("reader knows the following types:\n"); print(classify.ext(print.all=TRUE))
      return(NULL)
    }
  }
  if(typ==types[1])
  {
    ## .RData/.rda binary file
    file.out <- list()
    fns <- load(full.path)
    if(!quiet) { cat(paste(" found",paste(fns,collapse=","),"in",full.path,"\n")) }
    if (!is.null(want.type)) {
      for (ii in 1:length(fns)) {
        if(want.type %in% is(get(fns[ii])))
        { fns <- fns[ii] ; break }
      } }
    for (cc in 1:length(fns))
    {  file.out[[cc]] <- get(fns[cc]) }
    names(file.out) <- fns
    if(length(file.out)==1) { file.out <- file.out[[1]] }
  }
  if(typ==types[2])
  {
    # csv file
    file.out <- read.csv(full.path,stringsAsFactors=FALSE,...)
    file.out <- shift.rownames(file.out,override,warn=!quiet)
  }
  if(typ==types[3])
  {
    # other text .txt file
    detect <- suppressWarnings(get.delim(full.path,n=50,comment="#",large=10,one.byte=one.byte))
    if(length(detect)!=0) { def <- detect }
    first.10 <- readLines(full.path,n=10); hope10 <- length(first.10)
    if(hope10<3) { lown <- hope10 } else { lown <- 3 }
    if(hope10<2) { h2 <- lown <- 1 } else { h2 <- 2 }
    mini.parse <- strsplit(first.10,def)
    if(length(mini.parse)==0) { return(character(0)) }  # empty file
    splitto <- sapply(mini.parse,length)
    if(all(splitto[lown:hope10]==splitto[h2]) & (splitto[1]==(splitto[h2]-1))) {
      # was first row 1 delimiter shorter than all other rows? (i.e, header row)
      if(nchar(def)>1) {
        #read table only supports 1 byte delimiters (not regular expressions) 
        cat(" reading file with unspecified whitespace demiliting, check result imported correctly\n")
        cat(" header ignored!\n")
        ## START: SAME CHUNK AS BELOW FOR NON-HEADER VERSION ##
        raw.txt <- readLines(full.path); raw.txt <- raw.txt[raw.txt!=""]
        split.by.ws <- lapply(raw.txt[-1],strsplit,split=def) # hoping all same length
        lns <- sapply(split.by.ws,length)
        if(all(lns==1)) { split.by.ws <- Unlist(split.by.ws,1) } # if list got nested
        lns <- sapply(split.by.ws,length)
        if(length(unique(lns))!=1) {
          # some are different lengths, try removing blanks
          split.by.ws <- lapply(split.by.ws,function(X) { X[X!=""] })
          if(length(unique(sapply(split.by.ws,length)))!=1) {
            warning("attempt to read as white space implied different row lengths, returning list")
            return(split.by.ws)
          }
        }
        file.out <- t(as.df(split.by.ws))
        ## END: SAME CHUNK AS BELOW FOR NON-HEADER VERSION ##
        if(all(file.out[,1]=="") & ncol(file.out)>1) { file.out <- file.out[,-1] }
        rownames(file.out) <- NULL
      } else {
        file.out <- read.table(full.path, ..., sep=def,header=TRUE,row.names=1,stringsAsFactors=FALSE)
      }
    } else {
      ## probably has a header too!
      if(all(splitto==1)) {
        # probably just a vector file
        file.out <- readLines(full.path)
      } else {
        # some kind of delimited file
        if (all(splitto[1]==splitto[h2]))
        {
          # test whether there is a significant nchar difference between first and other rows
          # as a means of detecting whether the first row is a header row
          char.cnts <- sapply(mini.parse,nchar)
          #prv(char.cnts,mini.parse);
          lns <- sapply(char.cnts,length); 
          if(!sum(abs(diff(lns)))==0) {
            culprits <- lns!=Mode(lns)
            warning("Line(s): ",paste(which(culprits),collapse=","),"; seemed to have a different number of columns than the majority, will attempt to read anyway")
            mini.parse <- mini.parse[which(!culprits)]
            char.cnts <- sapply(mini.parse,nchar)
          }
          z.test <- function(X) { (X[1] - mean(X[-1],na.rm=TRUE))/max(1,sd(X[-1],na.rm=TRUE)) }
          Zs <- abs(apply(as.df(char.cnts),1,z.test))
          if(length(Zs)>0 & hope10>2) {
            #print(mini.parse); print(char.cnts); print(Zs)
            critZ <- abs(qnorm((h.test.p*(min(round(sqrt(length(Zs))),9)))/2))
            #cat("Zs mean",mean(Zs,na.rm=TRUE),"critZ",critZ,"\n")
            mZ <- mean(Zs,na.rm=TRUE)
           # prv(mZ,critZ,Zs,z.test,splitto,h.test.p,char.cnts,hope10)
            if(length(mZ)>0) {
             if(mZ<critZ) { hdr <- FALSE } else { hdr <- TRUE }
            } else { hdr <- TRUE }
            if(!is.na(header)) { if(is.logical(header)) { hdr <- header }}
            if(nchar(def)>1) {
              #read table only supports 1 byte delimiters (not regular expressions) 
              if(!quiet) {
                cat(" reading file with unspecified whitespace demiliting, check result imported correctly\n")
              }
              ## START: SAME CHUNK AS ABOVE FOR NON-HEADER VERSION ##
              raw.txt <- readLines(full.path); raw.txt <- raw.txt[raw.txt!=""]
              split.by.ws <- lapply(raw.txt,strsplit,split=def) # hoping all same length
              lns <- sapply(split.by.ws,length)
              if(all(lns==1)) { split.by.ws <- Unlist(split.by.ws,1) } # if list got nested
              lns <- sapply(split.by.ws,length)
              if(length(unique(lns))!=1) {
                # some are different lengths, try removing blanks
                split.by.ws <- lapply(split.by.ws,function(X) { X[X!=""] })
                if(length(unique(sapply(split.by.ws,length)))!=1) {
                  warning("attempt to read as white space implied different row lengths, returning list")
                  return(split.by.ws)
                }
              }
              file.out <- t(as.df(split.by.ws))
              ## END: SAME CHUNK AS ABOVE FOR NON-HEADER VERSION ##
              if(all(file.out[,1]=="") & ncol(file.out)>1) { file.out <- file.out[,-1] }
              if(nrow(file.out)>1 & hdr) {
                # assuming header row
                colnames(file.out) <- file.out[1,]; file.out <- file.out[-1,]
              }
              rownames(file.out) <- NULL
            } else {
              file.out <- read.delim(full.path, ..., sep=def,header=hdr,stringsAsFactors=FALSE)
              file.out <- shift.rownames(file.out,override,warn=!quiet)
            }
          } else {
            warning("file too small to determine structure, or other error in reader function, reverting to readLines")
            file.out <- readLines(full.path)
          }
        } else {
          if(!quiet) { 
            warning("*.txt file not delimited by",def)
            cat(" will just read as text\n")
          }
          file.out <- readLines(full.path)
        }
      }
    }
  }  
  if(exists("file.out")) {
    if(is.null(want.type)) {
      if(is.data.frame(file.out)) {
        if(is.null(rownames(file.out)) | (all(rownames(file.out)==paste(1:nrow(file.out)))) ) {
          if(dim(file.out)[2]==1) { 
            # if only 1 column, return a vector
            if(auto.vec) { file.out <- file.out[,1] }
          } else {
            # if first column has picked up default 1:nrow() rownames, remove it
            if(all(file.out[,1]==rownames(file.out))) { file.out[[1]] <- NULL }
          } } }
    }
  }
  return(file.out)
}


#' Classify file types readable by standard R I/O functions.
#' 
#' Look for known file extensions and classify as binary, comma-separated,
#' text format, or OTH=other; other files are assumed to be unreadable.
#' To read other files, need to specify more types manually.
#'
#' @param ext filenames or extensions to classify
#' @param more.txt more extensions that should be treated as txt
#' @param more.bin more extensions that should be treated as binary
#' @param more.csv more extensions that should be treated as csv
#' @param print.all setting to T, simply prints the list of supported ext
#' @return returns the 4 way classification for each file/extension
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{get.delim}}
#' @examples
#' classify.ext(c("test.txt","*.csv","tot","other","rda","test.RDatA"))
classify.ext <- function(ext=NULL,more.txt=NULL,more.bin=NULL,more.csv=NULL,print.all=FALSE) {
  bin <- c("RData","rda","bin")
  csv <- c("csv")
  txt <- c("txt","tab","map","vcf","bim","dat","fam","reg","cnv")
  if(is.character(more.txt)) { txt <- c(txt,more.txt) }
  if(is.character(more.bin)) { bin <- c(bin,more.bin) }
  if(is.character(more.csv)) { csv <- c(csv,more.csv) }
  bin <- tolower(bin); csv <- tolower(csv); Ext <- ext
  txt <- tolower(txt); ext <- tolower(ext);
  types <- c(bin,csv,txt)
  if(print.all) {  return(types) }
  if(length(grep(".",ext))>0) { 
    # maybe entered some filenames instead of extensions, try conversion
    ge <- get.ext(ext)
    ext[ge!=""] <- ge[ge!=""]
  }
  out <- rep("OTH",times=length(ext))
  out[ext %in% bin] <- "BIN"
  out[ext %in% txt] <- "TXT"
  out[ext %in% csv] <- "CSV"
  names(out) <- Ext
  return(out)
}


#' Shift the first column of a dataframe to rownames() if appropriate.
#' 
#' Checks whether the first column looks like IDs, and if so will.
#' remove the column, and move these values to rownames.
#'
#' @param dataf data.frame to run the conversion on
#' @param override assume col 1 is rownames, regardless of numeric() test
#' @param warn whether to display warnings if assumptions aren't met
#' @return returns vectors of strings of char, lengths X
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{reader}}
#' @examples
#' df1 <- data.frame(ID=paste("ID",101:110,sep=""),
#'                    scores=sample(70,10,TRUE)+30,age=sample(7,10,TRUE)+11)
#' shift.rownames(df1)
#' df2 <- data.frame(ID=paste(101:110),
#'                    scores=sample(70,10,TRUE)+30,age=sample(7,10,TRUE)+11)
#' shift.rownames(df2) # first col are all numbers, so no convert
#' shift.rownames(df2,override=TRUE) # override forces conversion
shift.rownames <- function(dataf,override=FALSE,warn=FALSE)
{
  if(is.data.frame(dataf))
  {
    if(ncol(dataf)<2) { return(dataf) } # only had 1 column, so can't convert it to rownames
    if(nrow(dataf)<1) { return(dataf) } # had no rows, so can't convert it to rownames
    typz <- sapply(sapply(dataf,is),"[[",1)
    rn <- paste(dataf[,1])
    dataout <- as.df(dataf[,-1])
    if (typz[1]=="character" & all(typz[-1] %in% c("numeric","integer")))
    {
      numerify <- TRUE
    } else {
      # test how many non-numeric in the first column (i.e. guess whether rownames)
      suppressWarnings(tst <- length(which(is.na(as.numeric(rn))))/nrow(dataout) )
      if (tst>=.75 | override)
      { 
        # test how many non-numeric in the new first column once rowname first col is removed
        suppressWarnings(tst <- length(which(is.na(as.numeric(paste(dataout[,1])))))/nrow(dataout) )
        if(tst<=.75) {
          numerify <- TRUE # convert rest of dataset to numeric (assume only 1st col was txt)
        } else {
          numerify <- FALSE }
      } else {
        if(warn) { warning("proposed rownames were mostly numbers") }
        return(dataf)
      }
    }    
    if(numerify) { for (dd in ncol(dataout))
    { suppressWarnings(dataout[,dd] <- as.numeric(as.character(dataout[,dd]))) }
    } 
    if(anyDuplicated(rn)) { if(warn) { warning("rownames not unique, so leaving as NULL") }; return (dataf) }
    rownames(dataout) <- rn
    return(dataout)
  } else {
    if(length(dim(dataf))==2)
    {
      if("character" %in% is(dataf[,1]))
      {
        if(length(unique(dataf[,1])==nrow(dataf)))
        {
          if(ncol(dataf)>1) {
            sup <- assess.dat.type(dataf)
            if(sup<2) { 
              if(warn){ warning("not sure if rownames should be added from col 1") }
              return(dataf)
            } else {
              rn <- dataf[,1]
              if(ncol(dataf)<3) { dataout <- as.df(dataf[,-1]) } else {
                dataout <- dataf[,-1] }
              if(anyDuplicated(rn)) { if(warn) { warning("rownames not unique, so leaving as NULL") }; return (dataf) }
              rownames(dataout) <- paste(rn)
              if(sup>=10) { for (dd in ncol(dataout))
              { suppressWarnings(dataout[,dd] <- as.numeric(as.character(dataout[,dd]))) }
              }    
              return(dataout)
            }            
          } else {
            if(warn) { warning("only a single column") }
            return(dataf)
          }
        } else {
          if(warn) { warning("duplicate row names") }
          return(dataf)
        }
      }
    } else {
      if(warn) { warning("must be formatted in rows and columns") }
      return(dataf)
    }
    if(warn) { warning("couldn't change first col, not a dataframe") }
    return(dataf)
  }
}


#' Internal function used by shift.rownames()
#' try to work out what format a file is in, whether IDs in column 1
assess.dat.type <- function(dat)
{
  support <- 0
  nchar.mn <- function(vec) { mean(nchar(paste(vec)),na.rm=TRUE) }
  if(length(unique(dat[,1])==nrow(dat)))
  { support <- support+1 }
  suppressWarnings(tst <- length(which(is.na(as.numeric(paste(dat[,1])))))/nrow(dat) )
  if (tst>=.75)
  { 
    suppressWarnings(tst <- length(which(is.na(as.numeric(paste(dat[,2])))))/nrow(dat) )
    if(tst<=.75) {
      support <- support+10
    } else {
      support <- support+1 }
  }
  if (ncol(dat)>2 & support <2) {
    # slow and potentially inaccurate so only run if in doubt
    rc <- dim(dat)
    if (nrow(dat)>100) { r.sel <- 1:100 } else { r.sel <- 1:rc[1] } 
    if (ncol(dat)>50) { c.sel <- 1:50 } else { c.sel <- 1:rc[2] } 
    col.chr <- apply(dat[r.sel,c.sel],2,nchar.mn)
    dif.to.1 <- mean(abs(col.chr[-1]-col.chr[1]))
    dif.to.0 <- mean(abs(col.chr[-1]-rev(col.chr[-1])))
    if (dif.to.0 < dif.to.1)
    { support <- support+1 } else { support <- support - 1 }
  }
  return(support)
}




#' Get the file extension from a file-name.
#'
#' @param fn filename(s) (with full path is ok too)
#' @return returns the (usually) 3 character file extension of a filename
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{rmv.ext}}
#' @examples
#' get.ext("/documents/nick/mydoc.xlsx")
#' get.ext(c("temp.cnv","temp.txt"))
get.ext <- function(fn) {
  # get file extension from a filename character string
  if(length(fn)<1) { warning("fn had length of zero"); return(fn) }
  if(all(is.na(fn)) | !is.character(fn)) { stop("fn should not be NA and should be of type character()") }
  strip.file.frags <- function(X) {
    file.segs <- strsplit(X,".",fixed=TRUE)[[1]]
    lss <- length(file.segs)
    if (lss>1) { out <- paste(file.segs[lss]) } else { out <- "" }
    return(out)
  }
  return(sapply(fn,strip.file.frags))
}

#' Remove the file extension from a file-name.
#'
#' Default is to only remove from a known list of file types,
#' this is to protect files with '.' which may not have an extension
#' This option can be changed, and more types can be specified too.
#'
#' @param fn filename(s) (with full path is ok too)
#' @param only.known logical, only remove extension if in the 'known' list
#' @param more.known character vector, add to the list of known extensions
#' @param print.known return the list of 'known' file extensions
#' @return returns the file name/path without the file extension
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{get.ext}}
#' @examples
#' rmv.ext(print.known=TRUE)
#' rmv.ext("/documents/nick/mydoc.xlsx")
#' rmv.ext(c("temp.cnv","temp.txt","temp.epi"))
#' # remove anything that looks like an extension
#' rmv.ext(c("temp.cnv","temp.txt","temp.epi"),only.known=FALSE) 
#' # add to list of known extensions
#' rmv.ext(c("temp.cnv","temp.txt","temp.epi"),more.known="epi") 
rmv.ext <- function(fn=NULL,only.known=TRUE,more.known=NULL,print.known=FALSE) {
  # remove file extension from a filename character string
  #if updating this function, also update internal copy in NCmisc
  known.ext <- c("TXT","RDATA","TAB","DAT","CSV","VCF","GCM","BIM","MAP","FAM",
                 "PFB","SH","R","CPP","H","DOC","DOCX","XLS","XLSX","PDF","JPG",
                 "BMP","PNG","TAR","GZ","CNV","PL","PY","ZIP","ORG","RDA","DSC","BCK",
                 "ABW","HTM","HTML",toupper(more.known))
  if(is.null(fn)) { 
    if(print.known) {
      return(known.ext)
    } else {
      warning("couldn't remove extension, not a character()"); return(fn) 
    }
  } else {
    if (all(is.na(fn))) { warning("couldn't remove extension, all values were NA"); return(fn) }
  }
  if(print.known) { cat("known file extensions:\n"); print(known.ext) }
  if(!is.character(fn)) { warning("couldn't remove extension, not a character()"); return(fn) }
  rmv.one <- function(X,known.ext) {
    file.segs <- strsplit(paste(X),".",fixed=TRUE)[[1]]
    lss <- length(file.segs)
    if (lss>1) { 
      if(only.known){
        if(toupper(file.segs[lss]) %in% known.ext) {
          out <- paste(file.segs[-lss],collapse=".") 
        } else { 
          out <- X
        }
      } else {
        out <- paste(file.segs[-lss],collapse=".") 
      }
    } else {
      out <- X 
    }
  }
  return(sapply(fn,rmv.one,known.ext=known.ext))
}


#' Internal function to assess whether data is a character or list of characters
is.ch <- function(x) { 
  # is function for character() or list of characters
  if(is.null(x)) { return(FALSE) }
  pt1 <- is.character(x)
  if(!pt1 & is.list(x)) { pt2 <- all(sapply(x,is.ch)) } else { pt2 <- pt1 }
  return(as.logical(pt1 | pt2))
}

#' Internal function 
as.df <- function(...) {
  ## unless 'stringsAsFactors' is called explicitly, wrap
  # as.data.frame so that stringsAsFactors is always FALSE
  test <- list(...)
  if(length(names(test))<1) { doit <- FALSE } else {
    if(names(test) %in% "stringsAsFactors") { doit <- TRUE } else { doit <- FALSE }
  }
  if(doit) {
    return(as.data.frame(...))
  } else {
    return(as.data.frame(...,stringsAsFactors=FALSE))
  }
}

#' Simple and robust way to create full-path file names.
#' 
#' Create a path with a file name, plus optional directory, prefix,
#' suffix, and file extension. dir/ext are robust, so that if they
#' already exist, the path produced will still make sense. Prefix is
#' applied after the directory, and suffix before the file extension.
#'
#' @param dir directory for the full path, if 'fn' already has a dir,
#'  then dir will be overridden. Auto add file separator if not present
#' @param fn compulsory vector of file names/paths
#' @param pref prefix to add in front of the file name
#' @param suf suffix to add after the file name, before the extension
#' @param ext file extension, will override an existing extension
#' @param must.exist the specified file must already exist, else error
#' @return returns vector of file names with the full paths
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' mydir <- "/Documents"
#' cat.path(mydir,"temp.doc")
#' # dir not added if one already present
#' cat.path(mydir,"/Downloads/me/temp.doc")
#' # using prefix and suffix
#' cat.path(mydir,"temp.doc","NEW",suf=5)
#' # changing the extension from .docx to .doc
#' cat.path(mydir,"temp.docx",ext="doc")
cat.path <- function(dir="",fn,pref="",suf="",ext="",must.exist=FALSE) 
{
  #if updating this function, also update internal copy in NCmisc
  dir.ch <- .Platform$file.sep
  if(is.list(fn) & is.ch(fn)) { fn <- unlist(fn) } #; 
  if(length(dir)>1) { dir <- dir[1]; cat("only first dir was used\n") }
  if(length(ext)>1) { ext <- ext[1]; cat("only first extension was used\n") }
  if(length(grep(dir.ch,fn))>0) {
    dir <- dirname(fn)  #split into dir and fn if fn has /'s
    fn <- basename(fn)
  }
  dir <- dir.force.slash(dir)
  if(ext!="") {
    #make sure ext includes the dot
    if(substr(ext,1,1)!=".")   { ext <- paste(".",ext,sep="") }
    #if ext is already built into suffix or filename, remove it from there
    fn <- rmv.ext(paste(fn))
    suf <- rmv.ext(paste(suf))
  }
  location <- paste(dir,pref,fn,suf,ext,sep="")
  if(any(!file.exists(location)) & must.exist) {
    warn <- paste("required file",location,"not found!")
    stop(warn)
  }
  return(location)
}



#' Internal function used by cat.path
dir.force.slash <- function(dir) {
  # make sure 'dir' directory specification ends in a / character
  if(!is.null(dim(dir))) { stop("dir should be a vector") }
  dir <- paste(dir)
  dir.ch <- .Platform$file.sep
  the.test <- (dir!="" & substr(dir,nchar(dir),nchar(dir))!=dir.ch)
  dir[the.test] <- paste(dir[the.test],dir.ch,sep="")
  return(dir)
}



#' Function to collect arguments when running R from the command line
#' 
#' Allows parameter specification by A=..., B=... in the command line
#' e.g, R < myScript.R M=1 NAME=John X=10.5, using commandArgs()
#'
#' @param arg.list the result of a commandArgs() call, or else NULL to
#'  initiate this call within the function
#' @param coms list of valid commands to look for, not case sensitive
#' @param def list of default values for each parameter (in same order)
#' @param verbose logical, whether to print to the console which assignments are made and warning messages
#' @param list.out logical, whether to return output as a list or data.frame 
#' @return returns dataframe showing the resulting values [column 1, "value"] for each 'coms' (rownames); or, if
#'  list.out=TRUE, then returns a list with names corresponding to 'coms' and values equivalent to 'value' column of 
#'  the data.frame that would be returned if list.out=FALSE
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' parse.args(c("M=1","NAME=John","X=10.5"),coms=c("M","X","NAME"))
#' parse.args(c("N=1")) # invalid command entered, ignored with warning
#' temp.fn <- "tempScript1234.R"
#' # make a temporary R Script file to call using the command line
#' # not run # writeLines(c("require(reader)","parse.args(coms=c('M','X','NAME'))"),con=temp.fn)
#' bash.cmd <- "R --no-save < tempScript1234.R M=1 NAME=John X=10.5"
#' # run above command in the terminal, or using 'system' below:
#' # not run # arg <- system(bash.cmd)
#' # not run # unlink(temp.fn) # delete temporary file
parse.args <- function(arg.list=NULL,coms=c("X"),def=0, list.out=F, verbose=TRUE)
{
  # parse arguments entered running R from the command line
  # using NAME=VALUE
  if(is.null(arg.list)) { arg.list <- commandArgs() }
  if(length(coms)>1 & length(def)==1) { def <- rep(def,length(coms)) }
  coms.original.case <- coms
  coms <- toupper(coms)
  outframe <- data.frame(value=paste(def),stringsAsFactors=FALSE)
  rownames(outframe) <- coms
  assign.cmds <- grep("=",arg.list,fixed=TRUE)
  if(length(assign.cmds)>0)
  {
    vars.lst <- strsplit(arg.list[assign.cmds],"=",fixed=TRUE)
    vars <- sapply(vars.lst,"[",1)
    vals <- sapply(vars.lst,tail,1)
    vals <- paste(vals)
    #vals <- as.integer(vals)
    if(any(toupper(vars) %in% coms))
    {
      which.coms <- match(toupper(vars),coms)
      for (cc in 1:length(which.coms))
      {
        if(!is.na(vals[cc]) & !is.na(which.coms[cc]))
        {
          if(verbose) { cat(paste("set",coms[which.coms[cc]],"=",vals[cc]),"\n") }
          outframe[coms[which.coms[cc]],1] <- paste(vals[cc])
        } else {
          if(verbose) {
            cat(paste(" skipping invalid variable",vars[cc],"or invalid value",vals[cc],"\n"))
          }
        }
      }
    } else {
      warning("command line arguments entered but none are valid")
    }
  } else {
    outframe <- NULL
  } 
  rownames(outframe) <- coms.original.case
  if(list.out) {
    outframe <- as.list(as.data.frame(t(outframe),stringsAsFactors=FALSE))
  }
  return (outframe)
}


#' Find the number of rows (lines) in a file.
#' 
#' Returns the number of lines in a file, which in the case of a datafile
#' will often correspond to the number of rows, or rows+1. Can also
#' do this for all files in the directory. File equivalent of nrow()
#'
#' @param fn name of the file(s) to get the length of
#' @param dir optional path for fn location, or specify all files in dir
#' @param all.in.dir select whether to extract length for all files in dir
#' @return returns length of file (or all files)
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{file.ncol}}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' write.table(matrix(rnorm(100),nrow=10),"temp.txt",col.names=FALSE)
#' file.nrow("temp.txt")
#' # use with caution, will be slow if dir contains large files
#' # not run # file.nrow(all.in.dir=TRUE) 
#' unlink("temp.txt")
#' setwd(orig.dir) # reset working directory to original
file.nrow <- function(fn="",dir="",all.in.dir=FALSE) {
  zip <- FALSE
  if(all(fn=="")) { 
    all.in.dir <- TRUE; fn <- paste(dir,"*",sep="") 
  } else {
    if(length(fn)>1) { warning("only first file used, to do multiple use all.in.dir") }
    if(!is.file(fn[1],dir)) { stop("Error: file did not exist") }
  }
  if(!check.linux.install("wc")) { return(suppressWarnings(wc.windows(fn))) }
  if(tolower(get.ext(fn))!="gz") { cmd <- paste("wc -l ",fn,sep="") } else {
    if(!all.in.dir) { cmd <- paste("zcat ",fn[1]," | wc -l ",sep="") ; zip <- TRUE } }
  linez <- system(cmd,intern=TRUE)
  dir.linez <- unique(c(grep("Is a directory",linez),grep("total",linez)))
  if(length(dir.linez)>0) { linez <- linez[-dir.linez] }
  splits <- strsplit(linez,"\t| +")
  splits <- lapply(splits,function(x) { x <- x[x!=""] })
  lens <- as.numeric(sapply(splits,"[",1))
  nms <- sapply(splits,"[",2)
  if(zip & all(is.na(nms))) { nms <- (fn[1]) }
  names(lens) <- nms
  return(lens)
}


# internal alternative for wc -l for windows
wc.windows <- function(fn) {
  dat.file <- file(fn)
  open(con=dat.file,open="r")
  eof <- FALSE; cc <- 0
  while(!eof) {
    eof <- length(readLines(dat.file,n=1))==0
    cc <- cc + 1
  }
  close(con=dat.file)
  return(cc-1)
}


#' Find the number of columns (lines) in a file.
#' 
#' Returns the number of columns in a datafile. File equivalent of ncol()
#'
#' @param fn name of the file(s) to get the length of
#' @param reader try to read the entire file to get a result, else
#'  looks at the top few lines (ignoring comments)
#' @param del specify a delimiter (else this will be auto-detected)
#' @param comment a comment symbol to ignore lines in files
#' @param skip number of lines to skip at top of file before processing
#' @param force try to read the file regardless of whether it looks
#'  like an invalid file type. Only use when you know the files are valid
#' @param excl.rn exclude rownames from column count (essentially subtract 1)
#' @return returns number of columns in file(s). If no delimiter, then =1
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{file.nrow}}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' write.table(matrix(rnorm(100),nrow=10),"temp.txt",col.names=FALSE,row.names=FALSE)
#' file.ncol("temp.txt",excl.rn=TRUE)
#' unlink("temp.txt")
#' # find ncol for all files in current directory:
#' # [NB: use with caution, will be slow if dir contains large files]
#' # not run # lf <- list.files(); if(length(lf)==0) { print("no files in dir") }
#' # lf <- lf[classify.ext(lf)=="TXT"]
#' # not run (only works if length(lf)>0) # file.ncol(lf) 
#' setwd(orig.dir) # reset working directory to original
file.ncol <- function(fn,reader=FALSE,del=NULL,comment="#",skip=0,force=FALSE,excl.rn=FALSE) {
  ## ncol function but for a file
  use.reader <- get("reader",mode="logical")
  if(!is.ch(fn)) { stop("Invalid format entered for filenames") }
  if(is.list(fn) | length(fn)>1 ) { return(sapply(fn,file.ncol,
      reader=reader,del=del,comment=comment,skip=skip,force=force)) }
  if(!is.file(fn,"")) { warning("Error: file did not exist"); return(NULL) }
  if(classify.ext(fn)=="BIN") { use.reader <- TRUE }
  if(classify.ext(fn)=="OTH") { 
    warn <- (paste(fn,"looks like an unsupported file type, may crash") )
    if(!force) { stop(warn) } else { warning(warn); return(NA) }
  }
  if(use.reader){
    fl <- reader(fn)
    if(is.null(dim(fl))) {
      if(is.vector(fl)) {
        fl <- force.frame(fl)
      } else {
        return(NA)
      }
    }
    ncols <- ncol(fl); rm(fl)
  } else {
    if(is.null(del)) {
      # suppress in case it's a vector file, e.g, ncol will be 1
      del <- suppressWarnings(get.delim(fn,n=6,comment=comment,skip=skip))
    }
    test.bit <- n.readLines(fn=fn,n=5,comment=comment,skip=skip)
    lens <- sapply(strsplit(test.bit,del),length)
    ncols <- ceiling(mean(lens,na.rm=TRUE))
  }
  if(excl.rn) { ncols <- ncols-1 }
  return(ncols)
}



#' Read 'n' lines (ignoring comments and header) from a file.
#' 
#' Useful when you don't know the length/structure of a file
#' and want a useful sample to look at. Can skip ahead in the file too.
#' Copes well when there are less than 'n' lines in the file.
#'
#' @param fn name of the file(s) to get the length of
#' @param n number of valid lines to attempt to read
#'  looks at the top few lines (ignoring comments)
#' @param comment a comment symbol to ignore lines in files
#' @param skip number of lines to skip at top of file before processing
#' @param header whether to allow for, and skip, a header row
#' @return returns the first n lines of the file meeting the criteria,
#'  or if 'skip' implies lines beyond the length of the file, the 
#'  result,will be truncated - although in this case, the last 
#'  line will always be read.
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' dat <- matrix(sample(100),nrow=10)
#' write.table(dat,"temp.txt",col.names=FALSE,row.names=FALSE)
#' n.readLines("temp.txt",n=2,skip=2,header=FALSE)
#' dat[3:4,]
#' unlink("temp.txt")
#' setwd(orig.dir) # reset working directory to original
n.readLines <- function(fn,n,comment="#",skip=0,header=TRUE)
{
  # read at least 'n' lines of a file, skipping lines and ignoring any starting with comment
  if(!file.exists(fn)) { warning("file doesn't exist"); return(NULL) }
  if(!is.character(comment)) { warning("illegal comment char, reverting to #"); comment <- "#" }
  rl <- 0; cc <- 0 + {if(is.numeric(skip)) skip else 0 }
  while(rl<n) { 
    test.bit <- readLines(fn,n+cc)
    if(skip>0 & length(test.bit>1)) { test.bit <- test.bit[-(1:(min((length(test.bit)-1),skip)))] }
    cmnt <- which(substr(test.bit,1,1)==comment)
    rl <- n+cc-length(cmnt)
    cc <- cc + length(cmnt)
  }
  if(length(cmnt)>0) { test.bit <- test.bit[-cmnt] } 
  if(length(test.bit)>1 & header) { test.bit <- test.bit[-1] }
  return(test.bit)
}


#' Determine the delimiter for a text data file.
#' 
#' Reads the first few lines of data in a text file and attempts to
#' infer what delimiter is in use, based on the 'delims' argument
#' that would result in the most consistent number of columns in the
#' first 'n' lines of data. Searches preferentially for delimiters
#' implying between 2 and 'large' columns, then for >large, and lastly
#' for 1 column if nothing else gives a match.
#'
#' @param fn name of the file to parse
#' @param n the number of lines to read to make the inference
#' @param comment a comment symbol to ignore lines in files
#' @param skip number of lines to skip at top of file before processing
#' @param delims the set of delimiters to test for
#' @param one.byte only check for one-byte delimiters, [e.g, whitespace regular expr is >1 byte]
#' @param large search initially for delimiters that imply more than 1, 
#'  and less than this 'large' columns; if none in this range, look next
#'  at >large.
#' @return returns character of the most likely delimiter
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{reader}}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' df <- data.frame(ID=paste("ID",101:110,sep=""),
#'   scores=sample(70,10,TRUE)+30,age=sample(7,10,TRUE)+11)
#' # save data to various file formats
#' test.files <- c("temp.txt","temp2.txt","temp3.csv")
#' write.table(df,file=test.files[1],col.names=FALSE,row.names=FALSE,sep="|",quote=TRUE)
#' write.table(df,file=test.files[2],col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)
#' write.csv(df,file=test.files[3])
#' # report the delimiters
#' for (cc in 1:length(test.files)) { 
#'   cat("\n",test.files[cc],": ")
#'   print(get.delim(test.files[cc])) }
#' unlink(test.files)
#' setwd(orig.dir) # reset working dir to original
get.delim <- function(fn,n=10,comment="#",skip=0,
                        delims=c("\t","\t| +"," ",";",","),large=10,one.byte=TRUE)  
{
  # test top 'n' lines to determine what delimeter the file uses
  if(!file.exists(fn)) { stop(paste("cannot derive delimiter as file",fn,"was not found"))}
  test.bit <- n.readLines(fn=fn,n=n,comment=comment,skip=skip)
  #print(test.bit)
  num.del <- list()
  if(any(nchar(delims)>1) & one.byte) { delims <- delims[-which(nchar(delims)>1)] }
  for (cc in 1:length(delims)) {
    fff <- nchar(delims[[cc]])==1
    num.del[[cc]] <- sapply(strsplit(test.bit,delims[[cc]],fixed=fff),length)
  }
  #prv(num.del)
  if(all(unlist(num.del)==1)) { 
    warning("not a delimited file, probably a vector file")
    return(NA)
  }
  # are there some delimiters that produce consistent ncol between rows?
  need.0 <- sapply(num.del,function(X) { sum(diff(X)) })
  num.del <- sapply(num.del,"[",1)
  if(any(!need.0)) {
    #rng <- range(num.del)
    candidates <- which(num.del>1 & num.del<=large & !need.0)
    #print(candidates)
    if(length(candidates)>0) { out <- candidates[1] 
    } else {
      candidates <- which(num.del>large & !need.0)
      if(length(candidates)>0) { out <- candidates[1]
      } else {
        candidates <- which(num.del==1 & !need.0)
        if(length(candidates)>0) { out <- candidates[1]
        } else {
          warning("no delimiters tried were able to produce a valid file spec")
          out <- NULL
        }
      }
    }
  } else {
    warning("no delimiters tried were able to produce a valid file spec")
    out <- NULL
  }
  #print(delims)
  return(delims[out])
}




#' Search for a directory to add to the path so that a file exists.
#' 
#' Looks for a file named 'fn' in 'dir', and if not found there, 
#' broadens the search to the list or vector of directorys, 'dirs'.
#' Returns the full path of the first match that exists.
#'
#' @param fn name of the file to search for
#' @param dir the first directory to look in (expected location)
#' @param dirs vector/list, a set of directories to look in should
#'  the file not be found in 'dir'.
#' @return if the file is found, returns the full path of the file,
#'  else returns an empty string ""
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{is.file}}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' l.fn <- "temp.txt"
#' writeLines("test",con=l.fn)
#' find.file(l.fn)
#' find.file(l.fn,dir=getwd())
#' unlink(l.fn)
#' # not run # common.places <- ## <<add local folder here>> ##
#' # not run # d.fn <- cat.path(common.places[1],l.fn)
#' # write this example file to the first of the folders #
#' # not run # if(!file.exists(d.fn)) {  writeLines("test2",con=d.fn) }
#' # search the local folders for a
#' # a file named 'temp.txt'
#' # not run # find.file(l.fn,dir=getwd(),dirs=common.places)
#' # unlink(d.fn) # run only if test file produced
#' setwd(orig.dir) # reset working dir to original
find.file <- function(fn,dir="",dirs=NULL) { 
  if(!is.ch(fn)) { return("") }
  for (cc in 1:length(fn)) {
    if(is.character(fn)) {
      fn[cc] <- add.dir.if.not(fn[cc],dir,dirs,TRUE) # e.g, add 'dir$ano' if needed
    } else {
      fn[[cc]] <- find.file(fn[[cc]],dir,dirs)
    }
  }
  return(fn)
}



#' Test whether a file exists in a target directory, or alternative
#' list of directories.
#' 
#' Looks for a file named 'fn' in 'dir', and if not found there, 
#' broadens the search to the list or vector of directorys, 'dirs'.
#' Returns TRUE or FALSE as to whether the file exists.
#'
#' @param fn name of the file to search for
#' @param dir the first directory to look in (expected location)
#' @param dirs vector/list, a set of directories to look in should
#'  the file not be found in 'dir'.
#' @param combine if a list is given, test whether ALL files valid
#' @return logical vector of whether each file was found,  or  if
#'  combine is true, then a single value whether ALL valid or not.
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{find.file}}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' l.fn <- "temp.txt"
#' writeLines("test",con=l.fn)
#' some.local.files <- narm(list.files()[1:10])
#' print(some.local.files)
#' is.file(l.fn)
#' is.file(l.fn,dir=getwd())
#' is.file(some.local.files)
#' # add a non-valid file to the list to see what happens
#' is.file(c(some.local.files,"fakefile.unreal"))
#' is.file(c(some.local.files,"fakefile.unreal"),combine=FALSE)
#' unlink(l.fn)
#' setwd(orig.dir) # reset working dir to original
is.file <- function(fn,dir="",dirs=NULL,combine=TRUE) {
  # if the path or raw filename 'fn' can refer to a file in dir/dirs
  # return TRUE, else FALSE
  if(combine) { FUN <- sapply } else { FUN <- lapply }
  X <- find.file(fn,dir,dirs)
  if(is.list(X)) {
    out <- FUN(X,is.file,dir,dirs,combine=combine)
  } else {
    if(is.character(X)) {
      out <- (X!="")
    } else {
      out <- rep(FALSE,times=length(X))  # was 'fn'??
    }
  }
  if(combine) {
    return(all(out))
  } else {
    return(out)
  }  
}


#' Internal function used by find.file
add.dir.if.not <- function(locs,dir="",dirs="",blank.if.not=TRUE,warn=FALSE) {
  # for any number of file names or paths, adds a directory if that is needed
  # to point to an existing file. can cycle through multiple directory possibilities
  # if desired; look first in 'dir', then try dirs if that fails
  if(is.list(dir)) { dir <- paste(unlist(dir)) }
  if(is.list(dirs)) { dirs <- paste(unlist(dirs)) }
  dir <- c(dir,dirs,"") # this makes sure 'dir' directory(s) are checked first, others, then current as last resort
  if(is.null(locs)) { 
    locs <- ""
    warning("filename was blank, might return nonsense location")
  } else {
    if(all(is.na(locs))) {
      locs <- ""
      warning("filename was NA, might return nonsense location")
    }
  }
  if(locs=="") { return("") } # otherwise will return first dir in list
  if(is.null(dir)) { dir <- "" ; if(warn) { warning("directory was blank") } }
  if(!is.character(locs)) { stop("Error: expecting filenames") }
  if(!is.character(dir)) { stop("Error: expecting directory names") }
  found <- FALSE
  for (cc in 1:length(locs)) {
    #if(!file.exists(locs[cc])) {
    dd <- 1
    while(!found & dd<=length(dir)) {
      test.next.dir <- cat.path(dir[dd],locs[cc])
      if(file.exists(test.next.dir))
      { 
        locs[cc] <- test.next.dir; found <- TRUE
      } else {
        found <- FALSE; dd <- dd + 1
      }
    }
    #}
    if(!file.exists(locs[cc]) & blank.if.not) { locs[cc] <- "" }
  }
  return(locs)
}




#' Convert a matrix or dataframe to fixed-width for nice file output
#' 
#' Pads each column to a common size so write.table() produces a
#' fixed width format that looks nice
#'
#' @param dat data.frame or matrix
#' @return returns dat with space padding as character
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk} #' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk} 
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' df <- data.frame(ID=paste("ID",99:108,sep=""),
#'   scores=sample(150,10,TRUE)+30,age=sample(16,10,TRUE))
#' dff <- conv.fixed.width(df)
#' write.table(df,file="notFW.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#' write.table(dff,file="isFW.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#' cat("Fixed-width:\n",paste(readLines("isFW.txt"),"\n"),sep="")
#' cat("standard-format:\n",paste(readLines("notFW.txt"),"\n"),sep="")
#' unlink(c("isFW.txt","notFW.txt"))
#' setwd(orig.dir) # reset working dir to original
conv.fixed.width <- function(dat) {
  # for a dataframe convert to fixed width format prior to writing to file
  padw <- function(X,L) { paste(spc(L-nchar(paste(X))),X,sep="") }
  if(is.null(dim(dat))) { stop("Error: need a matrix/dataframe as input") }
  for (cc in 1:ncol(dat)) {
    max.ch <- max(nchar(paste(dat[,cc])))
    dat[,cc] <- padw(dat[,cc],max.ch)
  }
  return(dat)
}

