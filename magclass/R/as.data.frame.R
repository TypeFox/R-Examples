
setMethod("as.data.frame",
  signature(x="magpie"),
  function(x,rev=1) 
  {
    if(rev==1) {
      dimnames(x)[[2]] <- getYears(x,as.integer=TRUE)
      if(is.null(dimnames(x)[[2]])) dimnames(x)[[2]] <- 0
      if(is.null(dimnames(x)[[3]])) dimnames(x)[[3]] <- "NA"
      if(any(dim(x)==0)) {
        return(data.frame())
      } else {
        x <- as.data.frame(as.table(x))
      }
      if(all(grepl(".",x[[3]],fixed=TRUE))) {
        tmp <- data.frame(t(matrix(unlist(strsplit(as.character(x[[3]]),split="\\.")),ncol=nrow(x))),stringsAsFactors=FALSE)
        for(i in 1:ncol(tmp)) tmp[[i]] <- factor(tmp[[i]],unique(tmp[[i]]))
        x <- cbind(x[,1:2],tmp,x[4])
      }
      colnames(x) <- c("Region","Year",paste("Data",1:(dim(x)[2]-3),sep=""),"Value")
      x <- cbind(Cell=suppressWarnings(as.integer(gsub("^[^\\.]*\\.","",x$Region))),x)
      x$Region <- gsub("\\..*$","",x$Region)
      return(x)
    } else if(rev==2) {
      x <- clean_magpie(x,what="sets")
      dimnames(x)[[2]] <- getYears(x,as.integer=TRUE)
      if(any(dim(x)==0)) {
        return(data.frame())
      } else {
        x <- as.data.frame(as.table(x), stringsAsFactors=FALSE)
      }
      names(x)[4] <- ".value"
      what <- ".value"
      types <- c(".spat",".temp",".data")
      for(i in 3:1) {
        if(grepl(".",names(x)[i],fixed=TRUE)) {
          tmp <- data.frame(t(matrix(unlist(strsplit(as.character(x[[i]]),split="\\.")),ncol=nrow(x))),stringsAsFactors=FALSE)
          names(tmp) <-  strsplit(names(x)[i],split="\\.")[[1]]
          if(i==1) {
            x <- cbind(tmp,x[2:length(x)])
          } else {
            x <- cbind(x[1:(i-1)],tmp,x[(i+1):length(x)])
          }
          what <- c(paste0(types[i],1:dim(tmp)[2]),what)
        } else {
          what <- c(paste0(types[i],1),what)
        }
      }
      #use other types than character if possible
      for(i in 1:ncol(x)) {
        if(is.character(x[[i]])) {
          x[[i]] <- type.convert(x[[i]],as.is=TRUE)
          if(is.character(x[[i]])) x[[i]] <- factor(x[[i]],unique(x[[i]]))
        }
      }
      attr(x,"dimtype") <- what
      return(x)
    } else {
      stop('Unknown revision "',rev,'"!')
    }
  }
)