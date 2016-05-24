#############################################################################################
## File: class.R
## Author: Xiaoyong Sun
## Date: 09/06/2010
## Goal: nonmem class
## Notes:
##      -
#############################################################################################


################################################################################
## Class
################################################################################
setClass("nonmem",
	representation(
          ## dir files
          file.cov="list",
          file.cor="list",
          file.coi="list",
          file.phi="list",
          
	        ## lst file
          file.lst="character",
          method="character",
          analysis="list",
          objt="character",
          objv="character",
          objs="character",
          
          ## tab file
          #file.tab="character",
          tabid="character",
          tabdata="data.frame"
          ),
	prototype=list(
          ## dir files
          file.cov=list(title=character(), data=data.frame()),
          file.cor=list(title=character(), data=data.frame()),
          file.coi=list(title=character(), data=data.frame()),
          file.phi=list(title=character(), data=data.frame()),
	
          ## lst file
          file.lst=character(),
          method=character(),
          analysis=list(),
          objt=character(),
          objv=character(),
          objs=character(),
          
          ## tab file
          #file.tab=character(),
          tabid=character(),
          tabdata=data.frame()
          )
)


setMethod("initialize", "nonmem", function(.Object,
          output.lst, output.tab, output.dir, delim=" ",...)
{
    ## check the first two arguments, the third one is optional
    if (missing(output.lst)) stop("\nNo outputfile!")
    if (!file.exists(output.lst)) stop("\nThere is no output.lst file in this path!")
    
    if (missing(output.tab)) stop("\nNo raw output file!")
    if (!file.exists(output.tab)) stop("\nThere is no output.tab file in this path!")
   
    ## pattern for standard output
    meth.pattern <- "#METH:"
    term.pattern <- "#TERM:"
    tere.pattern <- "#TERE:"
    objt.pattern <- "#OBJT:"
    objv.pattern <- "#OBJV:"
    objs.pattern <- "#OBJS:"
  
    ## pattern for additional output files
    cov.pattern <- ".cov"
    cor.pattern <- ".cor"
    coi.pattern <- ".coi"
    phi.pattern <- ".phi"
    
########################output lst##################################

    output <- readLines(output.lst, n=-1)
    file.standard <- output
    
    #METH
    meth.index <- grep(meth.pattern, output)
    
    #TERM, #TERE
    term.index <- grep(term.pattern, output)
    tere.index <- grep(tere.pattern, output)
    
    #OBJt
    objt.index <- grep(objt.pattern, output)
    
    #OBJV
    objv.index <- grep(objv.pattern, output)
    
    #OBJS
    objs.index <- grep(objs.pattern, output)

    if (length(meth.index) == 0) stop("\nNo method is found in output file!")

    meth <- NULL
    ter <- list()
    objt <- NULL
    objv <- NULL
    objs <- NULL
    for (i in 1:length(meth.index))
    {
        meth[i] <- strsplit(output[meth.index[i]], meth.pattern)[[1]][2]
        ter[[i]] <- output[(term.index[i]+1):(tere.index[i]-1)]
        objt[i] <- strsplit(output[objt.index[i]], objt.pattern)[[1]][2]

        #? maybe get from new raw output files
        objv[i] <- strsplit(output[objv.index[i]], objv.pattern)[[1]][2]
        objs[i] <- strsplit(output[objs.index[i]], objs.pattern)[[1]][2]
    }
    .Object@file.lst <- output
    .Object@method <- meth
    .Object@analysis <- ter
    .Object@objt <- objt
    .Object@objv <- objv
    .Object@objs <- objs
    
######################## output tab ##################################
     
    output <- readLines(output.tab,n=1)
    .Object@tabid <- output[1]
    tabdata <- read.table(output.tab, header=T, skip=1)
    .Object@tabdata <- tabdata

########################output dir##################################
    if (!missing(output.dir))
    {
        if (!file.exists(output.dir)) stop("\nThere is no output.dir folder in this path!")    
        all.filenames <- dir(path=output.dir)

        # cov
        cov.index <- grep(cov.pattern, all.filenames)
        cor.index <- grep(cor.pattern, all.filenames)
        coi.index <- grep(coi.pattern, all.filenames)
        phi.index <- grep(phi.pattern, all.filenames)
        
        covCheck <- length(cov.index) + length(cor.index) + length(coi.index) + length(phi.index)
        if (covCheck == 0) 
        {
          stop("\nThere is no .cov, .cor, .coi or .phi file in the output.dir folder. This option only works for NONMEM 7!")
        }
        
        if(length(cov.index) != 0)
        {
            filename <- paste(output.dir, all.filenames[cov.index], sep="/")
            cov.title <- readLines(filename, n=1)
            cov.data <- read.table(filename, header=T, skip=1)
            .Object@file.cov <- list(title=cov.title, data=cov.data)
        }

        if(length(cor.index) != 0)
        {
            filename <- paste(output.dir, all.filenames[cor.index], sep="/")
            cor.title <- readLines(filename, n=1)
            cor.data <- read.table(filename, header=T, skip=1)
            .Object@file.cor <- list(title=cor.title, data=cor.data)
        }

        if(length(coi.index) != 0)
        {
            filename <- paste(output.dir, all.filenames[coi.index], sep="/")
            coi.title <- readLines(filename, n=1)
            coi.data <- read.table(filename, header=T, skip=1)
            .Object@file.coi <- list(title=coi.title, data=coi.data)
        }

        if(length(phi.index) != 0)
        {
            filename <- paste(output.dir, all.filenames[phi.index], sep="/")
            phi.title <- readLines(filename, n=1)
            phi.data <- read.table(filename, header=T, skip=1)
            .Object@file.phi <- list(title=phi.title, data=phi.data)
        }
     }
    callNextMethod(.Object, ...)


})

setValidity("nonmem", function(object) 
{
     TRUE
})

if (is.null(getGeneric("non.lst"))) setGeneric("non.lst", function(object) standardGeneric("non.lst"))
if (is.null(getGeneric("non.lst.meth"))) setGeneric("non.lst.meth", function(object) standardGeneric("non.lst.meth"))
if (is.null(getGeneric("non.lst.term"))) setGeneric("non.lst.term", function(object) standardGeneric("non.lst.term"))
if (is.null(getGeneric("non.lst.objt"))) setGeneric("non.lst.objt", function(object) standardGeneric("non.lst.objt"))
if (is.null(getGeneric("non.lst.objv"))) setGeneric("non.lst.objv", function(object) standardGeneric("non.lst.objv"))
if (is.null(getGeneric("non.lst.objs"))) setGeneric("non.lst.objs", function(object) standardGeneric("non.lst.objs"))
if (is.null(getGeneric("non.tab"))) setGeneric("non.tab", function(object) standardGeneric("non.tab"))
if (is.null(getGeneric("non.cov"))) setGeneric("non.cov", function(object) standardGeneric("non.cov"))
if (is.null(getGeneric("non.cor"))) setGeneric("non.cor", function(object) standardGeneric("non.cor"))
if (is.null(getGeneric("non.coi"))) setGeneric("non.coi", function(object) standardGeneric("non.coi"))
if (is.null(getGeneric("non.phi"))) setGeneric("non.phi", function(object) standardGeneric("non.phi"))
if (is.null(getGeneric("non.select"))) setGeneric("non.select", function(object, ...) standardGeneric("non.select"))

setMethod("non.lst",signature(object="nonmem"), function(object) object@file.lst)
setMethod("non.lst.meth",signature(object="nonmem"), function(object) gsub("^ ", "", object@method))
setMethod("non.lst.term",signature(object="nonmem"), function(object) object@analysis)
setMethod("non.lst.objt",signature(object="nonmem"), function(object) object@objt) 
setMethod("non.lst.objv",signature(object="nonmem"), function(object) 
         {
              if (is.na(object@objv)) return(object@objv)         
              t1 <- unlist(strsplit(object@objv, " "))
              return(t1[t1!=""][2])
         })
setMethod("non.lst.objs",signature(object="nonmem"), function(object) 
         {
              if (is.na(object@objs)) return(object@objs)         
              t1 <- unlist(strsplit(object@objs, " "))
              return(t1[t1!=""][2])
         })
setMethod("non.tab",signature(object="nonmem"), function(object) list(ID=object@tabid, data=object@tabdata))

setMethod("non.cov",signature(object="nonmem"), function(object) object@file.cov)
setMethod("non.cor",signature(object="nonmem"), function(object) object@file.cor)
setMethod("non.coi",signature(object="nonmem"), function(object) object@file.coi)
setMethod("non.phi",signature(object="nonmem"), function(object) object@file.phi)

setMethod("non.select",signature(object="nonmem"), function(object, lines, sep=" ", ...) 
         {     
              if (missing(lines)) return(object$file.lst)         
              if (!is.numeric(lines)) stop("lines should be numeric!")
              if (length(lines)>0)
              {
                  t1 <- object@file.lst[lines]
                  t2 <- strsplit(t1, sep)
                  #deleteN <- function(a) a[a!=""]
                  t3 <- lapply(t2, function(a) a[a!=""])
                  tsum <- sum(diff(unlist(lapply(t3, function(a) length(a)))))
                  if (tsum==0) return( as.data.frame(do.call("rbind",t3)) )
                  else return(t3)
                  # as.numeric(as.character(f$V2))
                  
              }
              else
              {
                  stop("Length of lines is 0")
              }
              
         })




#non.f2s <- function(select.df)
#{
#    if( (!is.data.frame(select.df)) && (!is.matrix(select.df)) ) stop("The data must be data frame or matrix!")
#    options(scipen= -1000)
#    return(select.df)
#}

#stringAsFactor=FF...