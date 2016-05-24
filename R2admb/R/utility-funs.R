indent <- function(str,n=2) {
    paste(paste(rep(" ",n),collapse=""),str,sep="")
}

## format numbers equal width with leading zeros if necessary
numfmt <- function(x,len=length(x),sep=".") {
    paste(x,
          formatC(seq(len),width=format.info(seq(len)),flag="0"),
          sep=sep)
}

numfmt2 <- function(x,xdim,sep=".",sep2=".") {
    ff1 <- format.info(seq(xdim[1]))
    ff2 <- format.info(seq(xdim[2]))
    paste(x,
          outer(seq(xdim[1]),seq(xdim[2]),
                function(i,j) paste(formatC(i,width=format.info(ff1),flag="0"),
                                    formatC(j,width=format.info(ff2),flag="0"),sep=sep2)),
          sep=sep)
}

rep_pars <- function(parnames) {
    parnames <- as.character(parnames)
    parnames <- unlist(lapply(split(parnames,factor(parnames,levels=unique(parnames))),
                              function(x) {
                                  if (length(x)==1) x else numfmt(x)
                              }))
    parnames
}


str_contains <- function(x,y) {
    length(grep(x,y)>1)
}

get_names <- function(pars,info) {
    unlist(sapply(pars,
                  function(p) {
                      if (p %in% info$inits$vname) {
                          i <- match(p,info$inits$vname)
                          tt <- info$inits$type[i]
                          if (str_contains("number$",tt)) {
                              p
                          } else if (str_contains("vector$",tt)) {
                              numfmt(p,as.numeric(info$inits$X2[i]))
                              ## FIXME: may fail if this value needs to be parsed?
                          } else stop("can't handle matrix names yet")
                      } else if (p %in% info$raneff$vname) {
                          i <- match(p,info$raneff$vname)
                          numfmt(p,as.numeric(info$raneff$X2[i]))
                      }
                  }))
}




## summary() method ...
##  save model file with object???


read_chunk <- function(fn,sep="^#",maxlines=1000) {
    end <- FALSE
    ans <- character(maxlines)
    i <- 1
    has_sep <- function(x) length(grep(sep,x))>0
    while (!end) {
        tmp <- readLines(fn,n=1)
        if (i>1 && has_sep(tmp)) {
            end <- TRUE
            pushBack(tmp,fn)
        } else if (length(tmp)==0) {
            end <- TRUE
        } else {
            ans[i] <- tmp
            i <- i+1
        }
    }
    ans[1:(i-1)]
}

read_hst <- function(fn) {
    fn <- paste(fn,"hst",sep=".")
    if (!file.exists(fn)) {
        warning("file ",fn," not found: returning NULL")
        return(NULL)
    }
    f <- file(fn,open="r")
    r <- list()
    repeat {
        chunk <- read_chunk(f)
        ## cat(length(chunk),":",chunk[1],"\n") 
        if (length(chunk)<2) break
        r <- c(r,list(chunk))
    }
    labs <- sapply(r,"[[",1)
    ## single values
    w <- c(1:2,8,10)
    r[w] <- lapply(r[w],
                   function(x) as.numeric(x[2]))
    names(r)[w] <- c("sampsize","stepsize_scale","npars","rseed")
    ## vectors of value
    w <- c(3:7,9)
    r[w] <- lapply(r[w],function(x) 
        as.numeric(strsplit(gsub("^ +","",x[2])," ")[[1]]))
    names(r)[w] <- c("stepsizes","means","sdevs","lower","upper","mcmcparms")
    ## r$npars is NOT RELIABLE! use length(stepsizes instead)
    ## parameter matrices
    w <- 11:(10+length(r$stepsizes))
    r[w] <- lapply(r[w],
                   function(z) {
                       do.call(rbind,
                               lapply(z[-c(1,length(z))],
                                      function(x) {as.numeric(strsplit(x," ")[[1]])}))})
    names(r)[w] <- gsub("^#","",
                        gsub("\\[([0-9]+)\\]",".\\1",
                             gsub("; *//.*","",labs[w])))
    ans <- c(r[1:10],hists=list(r[w]))
    class(ans) <- "admb_hist"
    ans
}




## from Steve Martell
reptoRlist <- function(fn) {
    ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
    idx=sapply(as.double(ifile),is.na)
    vnam=ifile[idx] #list names
    nv=length(vnam) #number of objects
    A=list()
    r=0
    for(i in 1:nv) {
        ir=match(vnam[i],ifile)
        if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
        dum=NA
        if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
        if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrows=irr-ir-1,fill=TRUE))
        if(is.numeric(dum))#Logical test to ensure dealing with numbers
        {
            A[[ vnam[i ]]]=dum
        }
    }
    return(A)
}

strip_comments <- function(s) {
    ## strip comments (and terminal whitespace)
    gsub("[ \\\t]*//.*$","",s)
}

## processing variables
proc_var <- function(s,drop.first=TRUE,maxlen) {
    if (drop.first) s <- s[-1]
    ## drop LOCAL CALCS sections
    calclocs <- grep("_CALCS *$",s)
    if (length(calclocs)>0) {
        droplines <- unlist(apply(matrix(-calclocs,
                                         ncol=2,byrow=TRUE),1,function(x) seq(x[1],x[2])))
        s <- s[droplines]
    }
    ## strip comments & whitespace
    s2 <- gsub("^[ \\\t]*","",gsub("[;]*[ \\\t]*$","",strip_comments(s)))
    s2 <- s2[nchar(s2)>0]
    s2 <- s2[!grepl("+[ \\\t]*!!",s2)] ## strip !! lines
    words <- strsplit(s2," ")
    words <- lapply(words,function(x) x[x!=""])
    type <- sapply(words,"[[",1)
    rest <- sapply(words,"[[",2)
    rest2 <- strsplit(gsub("[(),]"," ",rest)," ")
    vname <- sapply(rest2,"[[",1)
    if (length(rest2) == 0) ret <- NULL else {
        maxlen0 <- max(sapply(rest2,length))
        if (missing(maxlen)) maxlen <- maxlen0
        else maxlen <- pmax(maxlen,maxlen0)
        opts <- t(sapply(rest2,
                         function(w) {
                             ## as.numeric()?
                             c(w[-1],rep(NA,maxlen+1-length(w)))
                         }))
        ret <- data.frame(type,vname,opts,stringsAsFactors=FALSE)
    }
    ret
}

drop_calcs <- function(s) {
    startcalc <- grep("^ *LOC(AL)*_CALCS",s)
    endcalc <- grep("^ *END_CALCS",s)
    ## calc may be ended by next section
    droplines <- numeric(0)
    for (i in seq_along(startcalc)) {
        if (length(endcalc)<i) endcalc[i] <- length(s)
        droplines <- c(droplines,startcalc[i]:endcalc[i])
    }
    if (length(droplines)>0) s <- s[-droplines]
    commcalc <- grep("^ +!!",s)
    if (length(commcalc)>0) s <- s[-commcalc]
    s
}

