setClass("magpie",contains="array",prototype=array(0,c(0,0,0)))

.dimextract <- function(x,i,dim,pmatch=FALSE,invert=FALSE) {
  if(length(i)==0) return(NULL)
  .countdots <- function(i) {
    return(nchar(gsub("[^\\.]","",i)))
  }
  if(.countdots(i[1])==.countdots(dimnames(x)[[dim]][1]) & pmatch==FALSE){
    #i vector seems to specify the full dimname
    if(!anyDuplicated(dimnames(x)[[dim]])) {
      if(invert) {
        return(which(!(dimnames(x)[[dim]] %in% i)))
      } else {
        match <- match(i,dimnames(x)[[dim]])
        if(any(is.na(match))) {
          stop("subscript out of bounds (\"",paste0(i[is.na(match)],collapse="\", \""),"\")")
        }
        return(match)
      }
    } else {
      warning("Your dimnames in dim=",dim," contain duplicates! This might lead to erronous results and bad code performance. Please try to avoid duplicates in dimnames under all circumstances!")
    }
  }
  
  pmatch1 <- ifelse(pmatch==TRUE | pmatch=="right",".*","")
  pmatch2 <- ifelse(pmatch==TRUE | pmatch=="left",".*","")
  tmp <- lapply(paste("(^|\\.)",pmatch1,.escapeRegex(i),pmatch2,"(\\.|$)",sep=""),grep,dimnames(x)[[dim]])
  if(any(vapply(tmp,length,length(tmp))==0)) stop("Data element(s) \"",paste(i[vapply(tmp,length,length(tmp))==0],collapse="\", \""),"\" not existent in MAgPIE object!")
  tmp <- unlist(tmp)
  if(invert) {
    tmp <- setdiff(1:dim(x)[dim],tmp)    
  }
  return(tmp)
}

.mselect_df <- function(x,df) {
  if(is.null(names(dimnames(x)))) stop("Dimnames must have names in order to use mselect!")
  dims <- dimCode(names(df),x)
  if(all(dims==0)) stop('None of the dimensions in the mapping could be found in the magpie object!')
  if(any(dims==0)) {
    dfmissing <- df[dims==0]
    df <- df[dims!=0]
    fdims <- dims
    dims <- dims[dims>0]
  } else {
    dfmissing <- NULL
  }
  if(anyDuplicated(dims)) stop('Dimension(s) "',paste(names(dims)[duplicated(dims)],collapse='", "'),'" appear(s) more than once in the given mapping!')
  
  if(any(dims<3)) {
    stop("Currently only mappings within the data dimensions are supported!")
  } else {
    sdims <- as.integer(round((dims-3)*10))
    maxdim <- nchar(gsub("[^\\.]","",names(dimnames(x))[3]))+1
    if(any(sdims>maxdim)) stop("Inconsistent dimension information. Data dimension specified which does not seem to exist!")
    if(nrow(df)>0) df <- matrix(sapply(df,.escapeRegex),dim(df),dimnames=dimnames(df))
    dmissing <- which(!(1:maxdim%in%sdims))
    sdims <- c(sdims,dmissing)
    for(d in dmissing) df <- cbind(df,"[^\\.]*")
    elems <- NULL
    search <-  paste0("^",apply(df[,sdims, drop=FALSE],1,paste,collapse="\\."),"$")
    found <- lapply(search,grep,getNames(x))
    x <- x[,,unlist(found)]
    length <- unlist(lapply(found,length))
    if(!is.null(dfmissing)) {
      if(length(dfmissing)>1) {
        name_extensions <- do.call("paste",c(dfmissing,sep="."))
      } else {
        name_extensions <- dfmissing[[1]]
      }
      getNames(x) <- paste(getNames(x),name_extensions[rep(1:length(name_extensions),length)],sep=".")
      getSets(x,fulldim=FALSE)[3] <- paste(getSets(x,fulldim=FALSE)[3],paste(names(dfmissing),collapse="."),sep=".")
    }
    if(any(length==0) & nrow(df)>0) {
      row_extensions <- gsub('\\.',".",sub('[^\\.]*','NA',sub("^\\^","",sub("\\$$","",search[length==0])),fixed=TRUE),fixed=TRUE)
      if(!is.null(dfmissing)) {
       row_extensions <- paste(row_extensions,name_extensions[length==0],sep=".") 
      }
      tmp <- new.magpie(getCells(x),getYears(x),row_extensions,0,sets=getSets(x))
      if(ndata(x)==0) {
        x <- tmp
      } else {
        x <- mbind(x,tmp)
      }
      if(getOption("magclass.verbosity")>1) cat("NOTE (.mselect_df): The following elements were added to x as they appeared in the mapping but not in x: ",paste0(row_extensions,collapse=", ")," (values set to 0)\n")
    }
    return(return(x))
  }
}

setMethod("[",
          signature(x = "magpie"),
          function (x, i, j, k, drop=FALSE,pmatch=FALSE,invert=FALSE) 
          {
            if(is.null(dim(x))) return(x@.Data[i])
            if(!missing(i)) {
              if(is.data.frame(i)) {
                return(.mselect_df(x,i))
              }
              if(is.factor(i)) i <- as.character(i)
              if(is.character(i)) i <- .dimextract(x,i,1,pmatch=pmatch,invert=invert)
            }
            if(!missing(j)) {
              if(is.factor(j)) j <- as.character(j)
              if(is.numeric(j) & any(j>dim(x)[2])) {
                j <- paste("y",j,sep="")
                if(invert) j <- getYears(x)[!(getYears(x) %in% j)]
              } else if(is.null(j)) {
                j <- 1:dim(x)[2]
              } else if(invert) {
                j <- getYears(x)[!(getYears(x) %in% j)]
              }
            }
            if(!missing(k)) {
              if(is.factor(k)) k <- as.character(k)
              if(is.character(k)) k <- .dimextract(x,k,3,pmatch=pmatch,invert=invert)
            }
            if(ifelse(missing(i),FALSE,is.array(i) | any(abs(i)>dim(x)[1]))) {
              #indices are supplied as array, return data as numeric
              return(x@.Data[i])
            } else if(missing(j) & ifelse(missing(k),TRUE,is.logical(k)) & ifelse(missing(i),FALSE,all(abs(i)<=dim(x)[1]))) {
              if(length(x@.Data[i,,,drop=FALSE])==0) {
                return(x@.Data[i])
              } else {
                x@.Data <- x@.Data[i,,,drop=FALSE]
                if(drop) x <- collapseNames(x)
                return(x)
              }
            } else {    
              if(!missing(k)) {
                if(is.logical(k)) {
                  # weird case in which k should be actually missing but gets the value of the next argument in the argument list (drop)
                  x@.Data <- x@.Data[i,j,,drop=FALSE] 
                } else {
                  x@.Data <- x@.Data[i,j,k,drop=FALSE]
                }
              } else {
                x@.Data <- x@.Data[i,j,,drop=FALSE]                
              }
              if(drop) x <- collapseNames(x)
              return(x)
            }
    }
)

setMethod("[<-",
          signature(x = "magpie"),
          function (x, i, j, k, value, pmatch=FALSE) 
          {       
            if(is.null(dim(x))) {
              tmp <- x@.Data
              tmp[i] <- k
              return(tmp)
            }
            if(!missing(i)) {
              if(is.factor(i)) i <- as.character(i)
              if(is.character(i)) i <- .dimextract(x,i,1,pmatch=pmatch) 
            }
            if(!missing(j)) {
              if(is.factor(j)) j <- as.character(j)
              if(is.numeric(j) & any(j>dim(x)[2])) j <- paste("y",j,sep="")
              else if(is.null(j)) j <- 1:dim(x)[2]
            }
            if(!missing(k)) {
              if(is.factor(k)) k <- as.character(k)
              if(is.character(k)) k <- .dimextract(x,k,3,pmatch=pmatch) 
            }
            if(missing(value)) {
              x@.Data[i] <- k 
              return(x)
            } else {
              if(is.magpie(value)){
                if(missing(i)) ii <- 1:dim(x)[1] else ii <- i
                if(missing(j)) jj <- 1:dim(x)[2] else jj <- j
                if(missing(k)) kk <- 1:dim(x)[3] else kk <- k
                value <- magpie_expand(value,x[ii,jj,kk])    
              } else if(length(value)!=length(x@.Data[i,j,k]) & length(value)!=1) {
                #dangerous writing of value as order might be wrong! 
                stop("Replacement does not work! Different replacement length!")
              } else if(length(value)!=1) {
                if(getOption("magclass.verbosity")>1) cat("NOTE ([<-): Dangerous replacement! As replacement value is not an MAgPIE object name checking is deactivated!\n")
              }
              x@.Data[i,j,k] <- value
              return(x)
            }
          }
)