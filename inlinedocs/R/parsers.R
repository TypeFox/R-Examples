do.not.generate <- structure(function
### Make a Parser Function used to indicate that certain Rd files
### should not be generated.
(...
### Character strings indicating Rd files without the .Rd suffix.
 ){
  filenames <- c(...)
  function(docs,...){
    for(fn in filenames){
      docs[[fn]] <- list()
    }
    docs$.overwrite <- TRUE
    docs
  }
### A Parser Function that will delete items from the outer
### Documentation List.
},ex=function(){
  silly.pkg <- system.file("silly",package="inlinedocs")
  owd <- setwd(tempdir())
  file.copy(silly.pkg,".",recursive=TRUE)

  ## define a custom Parser Function that will not generate some Rd
  ## files
  custom <- do.not.generate("silly-package","Silly-class")
  parsers <- c(default.parsers,list(exclude=custom))

  ## At first, no Rd files in the man subdirectory.
  man.dir <- file.path("silly","man")
  dir(man.dir)

  ## Running package.skeleton.dx will generate bare-bones files for
  ## those specified in do.not.generate, if they do not exist.
  package.skeleton.dx("silly",parsers)
  Rd.files <- c("silly-package.Rd","Silly-class.Rd","silly.example.Rd")
  Rd.paths <- file.path(man.dir,Rd.files)
  stopifnot(all(file.exists(Rd.paths)))
  
  ## Save the modification times of the Rd files
  old <- file.info(Rd.paths)$mtime

  ## make sure there is at least 2 seconds elapsed, which is the
  ## resolution for recording times on windows file systems.
  Sys.sleep(4) 
  
  ## However, it will NOT generate Rd for files specified in
  ## do.not.generate, if they DO exist already.
  package.skeleton.dx("silly",parsers)
  mtimes <- data.frame(old,new=file.info(Rd.paths)$mtime)
  rownames(mtimes) <- Rd.files
  mtimes$changed <- mtimes$old != mtimes$new
  print(mtimes)
  stopifnot(mtimes["silly-package.Rd","changed"]==FALSE)
  stopifnot(mtimes["Silly-class.Rd","changed"]==FALSE)
  stopifnot(mtimes["silly.example.Rd","changed"]==TRUE)

  unlink("silly",recursive=TRUE)
  setwd(owd)
})

### combine NULL objects.
combine.NULL<-function(x,y){
    if (class(x) == "NULL"){
        # print(paste("mm x=",x))
        # print(paste("mm class(x)=",class(x)))
        x=list("")
    }
    if (class(y) == "NULL"){
        # print(paste("mm y=",y))
        # print(paste("mm class(y)=",class(y)))
        y=list("")
    }
    return(combine(x,y))
}

### combine lists or character strings
combine <- function(x,y){
    UseMethod("combine")
}

### combine character strings by pasting them together
combine.character <- function(x,y)
    paste(x,y,sep="\n")

### combine lists by adding elements or adding to existing elements
combine.list <- function(x,y){
  toadd <- if(".overwrite"%in%names(y)){
    y <- y[names(y)!=".overwrite"]
    rep(TRUE,length(y))
  }else{
    !names(y)%in%names(x)
  }
  toup <- names(y)[!toadd]
  x[names(y)[toadd]] <- y[toadd]
  for(up in toup)x[[up]] <- combine(x[[up]],y[[up]])
  x
### A list, same type as x, but with added elements from y.
}


getSource <- function
### Extract a function's source code.
(fun.obj
### A function.
 ) {
      srcref <- attr(fun.obj, "srcref")
      if (!is.null(srcref)) {
        ##unlist(strsplit(as.character(srcref), "\n"))
        as.character(srcref)
      }
      else attr(fun.obj, "source")
### Source code lines as a character vector.
}

### Prefix for code comments used with grep and gsub.
prefix <- "^[ \t]*###[ \t]*"

decomment <- function
### Remove comment prefix and join lines of code to form a
### documentation string.
(comments
### Character vector of prefixed comment lines.
 ){
  gsub(prefix,"",comments)
### String without prefixes or newlines.
}

forall <- function
### For each object in the package that satisfies the criterion
### checked by subfun, parse source using FUN and return the resulting
### documentation list.
(FUN,
### Function to apply to each element in the package.
 subfun=function(x)TRUE
### Function to select subsets of elements of the package, such as
### is.function. subfun(x)==TRUE means FUN will be applied to x and
### the result will be returned.
 ){
  FUN <- FUN
  f <- function(objs,docs,...){
    if(length(objs)==0)return(list())
    objs <- objs[sapply(objs,subfun)]
    L <- list()
    on.exit(cat(sprintf("Parser Function failed on %s\n",N)))
    for(N in union(names(docs),names(objs))){
      o <- objs[[N]]
      L[[N]] <- FUN(src=getSource(o),
                    name=N,objs=objs,o=o,docs=docs,doc=docs[[N]],...)
    }
    on.exit()## remove warning message
    L
  }
  class(f) <- c("allfun","function")
  f
### A Parser Function.
}

### Print method for functions constructed using forall.
print.allfun <- function(x,...){
  e <- environment(x)
  cat("Function to apply to every element.\nselector:")
  print(e$subfun)
  cat("processor:")
  print(e$FUN)
}

### For each function in the package, do something.
forfun <- function(FUN)forall(FUN,is.function)

kill.prefix.whitespace <- function
### Figure out what the whitespace preceding the example code is, and
### then delete that from every line.
(ex
### character vector of example code lines.
 ){
  tlines <- gsub("\\s*","",ex)
  ##tlines <- gsub("#.*","",tlines)
  prefixes <- unique(gsub("\\S.*","",ex[tlines!=""]))
  FIND <- prefixes[which.min(nchar(prefixes))]
  ## Eliminate leading tabulations or 2/4 spaces
  sub(FIND, "", ex)
### Character vector of code lines with preceding whitespace removed.
}

prefixed.lines <- structure(function(src,...){
### The primary mechanism of inline documentation is via consecutive
### groups of lines matching the specified prefix regular expression
### "\code{^### }" (i.e. lines beginning with "\code{### }") are
### collected as follows into documentation sections:\describe{
### \item{description}{group starting at line 2 in the code}
### \item{arguments}{group following each function argument}
### \item{value}{group ending at the penultimate line of the code}}
### These may be added to by use of the \code{##<<} constructs
### described below.
  clines <- grep(prefix,src)
  if(length(clines)==0)return(list())
  bounds <- which(diff(clines)!=1)
  starts <- c(1,bounds+1)
  ends <- c(bounds,length(clines))
  ## detect body of function using paren matching
  code <- gsub("#.*","",src)
  f <- function(ch)cumsum(nchar(gsub(sprintf("[^%s]",ch),"",code)))
  parens <- f("(")-f(")")
  body.begin <- which(diff(parens)<0 & parens[-1]==0)+2
  if(length(body.begin)==0)body.begin <- 1 ## rare cases
  is.arg <- function(){
    gres <- grep("^\\s*#",src[start-1],perl=TRUE)
    0 == length(gres) && start<=body.begin
  }
  res <- list()
  for(i in seq_along(starts)){
    start <- clines[starts[i]]
    end <- clines[ends[i]]
    processed <- gsub("#.*","",gsub("[ }]","",src[(end+1):length(src)]))
    lab <- if(all(processed==""))"value"
    else if(start==2)"description"
    else if(is.arg()){
      ##twutz: strip leading white spaces and brackets and ,
      arg <- gsub("^[ \t(,]*", "", src[start - 1])
      arg <- gsub("^([^=,]*)[=,].*", "\\1", arg)
      ##twutz: remove trailing whitespaces
      arg <- gsub("^([^ \t]*)([ \t]+)$","\\1",arg)
      arg <- gsub("...", "\\dots", arg, fixed = TRUE)
      paste("item{",arg,"}",sep="")
    } else {
      next;
    }
    res[[lab]] <- decomment(src[start:end])
  }
  res
},ex=function(){
test <- function
### the description
(x,
### the first argument
 y ##<< another argument
 ){
  5
### the return value
##seealso<< foobar
}
src <- getSource(test)
prefixed.lines(src)
extract.xxx.chunks(src)
})

extract.xxx.chunks <- function # Extract documentation from a function
### Given source code of a function, return a list describing inline
### documentation in that source code.
(src,
### The source lines of the function to examine, as a character
### vector.
 name.fun="(unnamed function)",
### The name of the function/chunk to use in warning messages.
 ...
### ignored.
 ){
  res <- list()
  ##details<< For simple functions/arguments, the argument may also be
  ## documented by appending \code{##<<} comments on the same line as the
  ## argument name. Mixing this mechanism with \code{###} comment lines for
  ## the same argument is likely to lead to confusion, as the \code{###}
  ## lines are processed first.
  #arg.pat <- paste("^[^=,#]*?([\\w\\.]+)\\s*([=,].*|\\)\\s*)?",
  #                 "<<\\s*(\\S.*?)\\s*$",
  #                 sep="##") # paste avoids embedded trigger fooling the system
   #tw: removed first comma
   arg.pat <- paste("^[^=#]*?([\\w\\.]+)\\s*([=,].*|\\)\\s*)?",
	   "<<\\s*(\\S.*?)\\s*$",
   		sep="##") # paste avoids embedded trigger fooling the system

  skeleton.fields <- c("alias","details","keyword","references","author",
                       "note","seealso","value","title","description",
                       "describe","end")
  ##details<< Additionally, consecutive sections of \code{##} comment
  ## lines beginning with \code{##}\emph{xxx}\code{<<} (where
  ## \emph{xxx} is one of the fields: \code{alias}, \code{details},
  ## \code{keyword}, \code{references}, \code{author}, \code{note},
  ## \code{seealso}, \code{value}, \code{title} or \code{description})
  ## are accumulated and inserted in the relevant part of the .Rd
  ## file.
  ##
  ## For \code{value}, \code{title}, \code{description} and function
  ## arguments, these \emph{append} to any text from "prefix"
  ## (\code{^### }) comment lines, irrespective of the order in the
  ## source.
  ##
  ## When documenting S4 classes, documentation from \code{details}
  ## sections will appear under a section \code{Objects from the Class}. That
  ## section typically includes information about construction methods
  ## as well as other description of class objects (but note that the
  ## class Slots are documented in a separate section).

  ## but this should not appear, because separated by a blank line
  extra.regexp <- paste("^\\s*##(",paste(skeleton.fields,collapse="|"),
                        ")<<\\s*(.*)$",sep="")
  cont.re <- "^\\s*##\\s*"
  in.describe <- 0
  first.describe <- FALSE
  k <- 1
  in.chunk <- FALSE
  end.chunk <- function(field,payload)
    {
      if ( "alias" == field ){
        ##note<< \code{alias} extras are automatically split at new lines.
        payload <- gsub("\\n+","\\}\n\\\\alias\\{",payload,perl=TRUE)
        chunk.sep <- "}\n\\alias{"
      } else if ( "keyword" == field ){
        ##keyword<< documentation utilities
        ##note<< \code{keyword} extras are automatically split at white space,
        ## as all the valid keywords are single words.
        payload <- gsub("\\s+","\\}\n\\\\keyword\\{",payload,perl=TRUE)
        chunk.sep <- "}\n\\keyword{"
      } else if ( "title" == field ){
        chunk.sep <- " "
      } else if ( "description" == field ){
        chunk.sep <- "\n"
      } else {
        ##details<< Each separate extra section appears as a new
        ## paragraph except that: \itemize{\item empty sections (no
        ## matter how many lines) are ignored;\item \code{alias} and
        ## \code{keyword} sections have special rules;\item
        ## \code{description} should be brief, so all such sections
        ## are concatenated as one paragraph;\item \code{title} should
        ## be one line, so any extra \code{title} sections are
        ## concatenated as a single line with spaces separating the
        ## sections.}
        chunk.sep <- "\n\n"
      }
      chunk.res <- NULL
      if ( !grepl("^\\s*$",payload,perl=TRUE) )
        chunk.res <-
          if ( is.null(res[[field]]) ) payload
          else paste(res[[field]], payload, sep=chunk.sep)
      invisible(chunk.res)
    }
  while ( k <= length(src) ){
    line <- src[k]
    ##print(line)
    ##if(grepl("^$",line))browser()
    if ( grepl(extra.regexp,line,perl=TRUE) ){
      ## we have a new extra chunk - first get field name and any payload
      new.field <- gsub(extra.regexp,"\\1",line,perl=TRUE)
      new.contents <- gsub(extra.regexp,"\\2",line,perl=TRUE)
      ##cat(new.field,"\n-----\n",new.contents,"\n\n")
      ##details<< As a special case, the construct \code{##describe<<} causes
      ## similar processing to the main function arguments to be
      ## applied in order to construct a describe block within the
      ## documentation, for example to describe the members of a
      ## list. All subsequent "same line" \code{##<<} comments go into that
      ## block until terminated by a subsequent \code{##}\emph{xxx}\code{<<} line.
      if ( "describe" == new.field ){
        ##details<< Such regions may be nested, but not in such a way
        ## that the first element in a \code{describe} is another
        ## \code{describe}.  Thus there must be at least one
        ## \code{##<<} comment between each pair of
        ## \code{##describe<<} comments.
        if ( first.describe ){
          stop("consecutive ##describe<< at line",k,"in",name.fun)
        } else {
          if ( nzchar(new.contents) ){
            if ( is.null(payload) || 0 == nzchar(payload) ){
              payload <- new.contents
            } else {
              payload <- paste(payload,new.contents,sep="\n\n")
            }
          }
          first.describe <- TRUE
        }
      } else if ( "end" == new.field ){
        ##details<< When nested \code{describe} blocks are used, a comment-only
        ## line with \code{##end<<} terminates the current level only; any
        ## other valid \code{##}\emph{xxx}\code{<<} line terminates
        ## all open describe blocks.
        if ( in.describe>0 ){
          ## terminate current \item and \describe block only
          if ( "value" == cur.field && 1 == in.describe ){
            payload <- paste(payload,"}",sep="")
          } else {
            payload <- paste(payload,"}\n}",sep="")
          }
          in.describe <- in.describe-1;
        } else {
          warning("mismatched ##end<< at line ",k," in ",name.fun)
        }
        if ( nzchar(new.contents) ){
          if ( nzchar(payload) ){
            payload <- paste(payload,new.contents,sep="\n")
          } else {
            payload <- new.contents
          }
        }
      } else {
        ## terminate all open \describe blocks (+1 because of open item)
        if ( 0 < in.describe ){
          if ( "value" != cur.field ){  # value is implicit describe block
            payload <- paste(payload,"}",sep="")
          }
          while ( in.describe>0 ){
            payload <- paste(payload,"}",sep="\n")
            in.describe <- in.describe-1;
          }
        }
        ## finishing any existing payload
        if ( in.chunk ) res[[cur.field]] <- end.chunk(cur.field,payload)
        in.chunk <- TRUE
        cur.field <- new.field
        payload <- new.contents
        ##note<< The "value" section of a .Rd file is implicitly a describe
        ## block and \code{##}\code{value}\code{<<} acts accordingly. Therefore
        ## it automatically enables the describe block itemization (##<< after
        ## list entries).
        if ( "value" == new.field ){
          first.describe <- TRUE;
        }
      }
    } else if ( in.chunk && grepl(cont.re,line,perl=TRUE) ){
      ## append this line to current chunk
      if ( !grepl(prefix,line,perl=TRUE) ){
        ##describe<< Any lines with "\code{### }" at the left hand
        ## margin within the included chunks are handled separately,
        ## so if they appear in the documentation they will appear
        ## before the \code{##}\emph{xxx}\code{<}\code{<} chunks.
### This one should not appear.
        stripped <- gsub(cont.re,"",line,perl=TRUE)
        if ( nzchar(payload) ){
          payload <- paste(payload,stripped,sep="\n")
        } else {
          payload <- stripped
        }
      }
    } else if ( grepl(arg.pat,line,perl=TRUE) ){
      not.describe <- (0==in.describe && !first.describe)
      if ( in.chunk && not.describe){
        res[[cur.field]] <- end.chunk(cur.field,payload)
      }
      comment <- gsub(arg.pat,"\\3",line,perl=TRUE);
      arg <- gsub(arg.pat,"\\\\item\\{\\1\\}",line,perl=TRUE)
      in.chunk <- TRUE
      if ( not.describe ){
        ## TDH 2010-06-18 For item{}s in the documentation list names,
        ## we don't need to have a backslash before, so delete it.
        arg <- gsub("^[\\]+","",arg)
        cur.field <- gsub("...","\\dots",arg,fixed=TRUE) ##special case for dots
        payload <- comment
      } else {
        ## this is a describe block, so we need to paste with existing
        ## payload as a new \item.
        if ( first.describe ){
          ## for first item, need to add describe block starter
          if ( "value" == cur.field ){
            payload <- paste(payload,"\n",arg,"{",sep="")
          } else {
            payload <- paste(payload,"\\describe{\n",arg,"{",sep="")
          }
          first.describe <- FALSE
          in.describe <- in.describe+1
        } else {
          ## subsequent item - terminate existing and start new
          payload <- paste(payload,"}\n",arg,"{",sep="")
        }
        if ( nzchar(comment) ){
          payload <- paste(payload,comment,sep="")
        }
      }
    } else if ( in.chunk ){
      if ( 0 == in.describe && !first.describe ){
        ## reached an end to current field, but need to wait if in.describe
        res[[cur.field]] <- end.chunk(cur.field,payload)
        in.chunk <- FALSE
        cur.field <- NULL
        payload <- NULL
      }
    }
    k <- k+1
  }
  ## finishing any existing payload
  if ( 0 < in.describe ){
    if ( "value" != cur.field ){    # value is implicit describe block
      payload <- paste(payload,"}",sep="")
    }
    while ( in.describe>0 ){
      payload <- paste(payload,"}",sep="\n")
      in.describe <- in.describe-1;
    }
  }
  if ( in.chunk ) res[[cur.field]] <- end.chunk(cur.field,payload)
  res
### Named list of character strings extracted from comments. For each
### name N we will look for N\{...\} in the Rd file and replace it
### with the string in this list (implemented in modify.Rd.file).
}

leadingS3generic <- function # check whether function name is an S3 generic
### Determines whether a function name looks like an S3 generic function
(name,                     ##<< name of function
 env,                      ##<< environment to search for additional generics
 ...)                      ##<< ignored here
{
  ##details<< This function is one of the default parsers, but exposed as
  ## possibly of more general interest. Given a function name of the form
  ## x.y.z it looks for the generic function x applying to objects of class
  ## y.z and also for generic function x.y applying to objects of class z.
  ##
  parts <- strsplit(name, ".", fixed = TRUE)[[1]]
  l <- length(parts)
  if (l > 1) {
    for (i in 1:(l - 1)) {
      ## Look for a generic function (known by the system or defined
      ## in the package) that matches that part of the function name
      generic <- paste(parts[1:i], collapse = ".")
      if (any(generic %in% getKnownS3generics()) ||
          findGeneric(generic, env) != "") {
        object <- paste(parts[(i + 1):l], collapse = ".")
        ##details<< Assumes that the first name which matches any known
        ## generics is the target generic function, so if both x and x.y
        ## are generic functions, will assume generic x applying to objects
        ## of class y.z
        ##value<< If a matching generic found returns a list with a single component:
        return(list(.s3method=c(generic, object))) ##<< a character vector containing generic name and object name.
      }
    }
  }
  ##value<< If no matching generic functions are found, returns an empty list.
  list()
}

### Parsers for each function that are constructed automatically. This
### is a named list, and each element is a parser function for an
### individual object.
forfun.parsers <-
  list(prefixed.lines=prefixed.lines,
       extract.xxx.chunks=extract.xxx.chunks,
       ## title from first line of function def
       title.from.firstline=function(src,...){
         first <- src[1]
         if(!is.character(first))return(list())
         if(!grepl("#",first))return(list())
         list(title=gsub("[^#]*#\\s*(.*)","\\1",first,perl=TRUE))
       },
       ## PhG: it is tests/FUN.R!!! I would like more flexibility here
       ## please, let me choose which dir to use for examples!
       ## Get examples for FUN from the file tests/FUN.R
       examples.from.testfile=function(name,...){
         tsubdir <- getOption("inlinedocs.exdir")
         if (is.null(tsubdir)) tsubdir <- "tests"	# Default value
         tfile <- file.path("..",tsubdir,paste(name,".R",sep=""))
         if(file.exists(tfile))
           list(examples=readLines(tfile))
         else list()
       },
       definition.from.source=function(doc,src,...){
         def <- doc$definition
         is.empty <- function(x)is.null(x)||x==""
         if(is.empty(def) && !is.empty(src))
           list(definition=src)
         else list()
       })

### List of Parser Functions that can be applied to any object.
forall.parsers <-
  list(## Fill in author from DESCRIPTION and titles.
       author.from.description=function(desc,...){
         list(author=desc[,"Author"])
       },
       ## The format section sometimes causes problems, so erase it.
       erase.format=function(...){
         list(format="")
       },
       ## Convert the function name to a title.
       title.from.name=function(name,doc,...){
         if("title"%in%names(doc))list() else
         list(title=gsub("[._]"," ",name))
       },
       ## PhG: here is what I propose for examples code in the 'ex' attribute
       examples.in.attr =  function (name, o, ...) {
         ex <- attr(o, "ex")
         if (!is.null(ex)) {
           ## Special case for code contained in a function
           if (inherits(ex, "function")) {
             ## If source is available, start from there
             src <- getSource(ex)
             if (!is.null(src)) {
               ex <- src
             } else { ## Use the body of the function
               ex <- deparse(body(ex))
             }
             ## Eliminate leading and trailing code
             ex <- ex[-c(1, length(ex))]
             if(ex[1]=="{")ex <- ex[-1]
             ## all the prefixes
             ex <- kill.prefix.whitespace(ex)
             ## Add an empty line before and after example
             ex <- c("", ex, "")
           }
           list(examples = ex)
         } else list()
       },collapse=function(doc,...){
         L <- lapply(doc,paste,collapse="\n")
         L$.overwrite <- TRUE
         L
       },tag.s3methods=leadingS3generic
       )

### List of parser functions that operate on single objects. This list
### is useful for testing these functions.
lonely <- structure(c(forall.parsers,forfun.parsers),ex=function(){
  f <- function # title
### description
  (x, ##<< arg x
   y
### arg y
   ){
    ##value<< a list with elements
    list(x=x, ##<< original x value
         y=y, ##<< original y value
         sum=x+y) ##<< their sum
    ##end<<
  }
  src <- getSource(f)
  lonely$extract.xxx.chunks(src)
  lonely$prefixed.lines(src)
})

extra.code.docs <- function # Extract documentation from code chunks
### Parse R code to extract inline documentation from comments around
### each function. These are not able to be retreived simply by
### looking at the "source" attribute. This is a Parser Function that
### can be used in the parser list of package.skeleton.dx(). TODO:
### Modularize this into separate Parsers Functions for S4 classes,
### prefixes, ##<<blocks, etc. Right now it is not very clean!
(code,
### Code lines in a character vector containing multiple R objects to
### parse for documentation.
 objs,
### The objects defined in the code.
 ...
### ignored
 ){
  parsed <- extract.file.parse(code)
  extract.docs.try <- function(o,on)
    {
      ## Note: we could use parsed information here too, but that
      ## would produce different results for setMethodS3 etc.
      doc <- list()
      if ( !is.null(parsed[[on]]) ){
        if ( !is.na(parsed[[on]]@code[1]) ){ # no code given for generics
          doc$definition <- paste(parsed[[on]]@code)
        }
        if(!"description"%in%names(doc) && !is.na(parsed[[on]]@description) ){
          doc$description <- parsed[[on]]@description
        }
        ## if ( "setMethodS3" == parsed[[on]]@created ){
        ##   gen <- leadingS3generic(on,topenv())
        ##   if ( 0 < length(gen) ){
        ##     doc$.s3method <- gen$.s3method
        ##     cat("S3method(",gen$.s3method[1],",",gen$.s3method[2],")\n",sep="")
        ##   }
        ## }
      }
      if("title" %in% names(doc) && !"description" %in% names(doc) ){
        ## For short functions having both would duplicate, but a
        ## description is required. Therefore automatically copy title
        ## across to avoid errors at package build time.
        doc$description <- doc$title
      }
      doc
    }
  extract.docs <- function(on){
    res <- try({o <- objs[[on]]
                extract.docs.try(o, on)},FALSE)
    if(class(res)=="try-error"){
      cat("Failed to extract docs for: ",on,"\n\n")
      list()
    } else if(0 == length(res) && inherits(objs[[on]],"standardGeneric")){
      NULL
    } else if(0 == length(res) && "function" %in% class(o)
              && 1 == length(osource <- getSource(o))
              && grepl(paste("UseMethod(",on,")",sep="\""),osource)
              ){
      ## phew - this should only pick up R.oo S3 generic definitions like:
      ## attr(*, "source")= chr "function(...) UseMethod(\"select\")"
      NULL
    } else res
  }
  doc.names <- names(objs)
  res <- sapply(doc.names,extract.docs,simplify=FALSE)
  ## Special processing for S4 classes as they do not appear in normal ls()
  for ( nn in names(parsed) ){
    if ( parsed[[nn]]@created == "setClass" ){
      S4class.docs <- extract.docs.setClass(parsed[[nn]])
      docname <- paste(nn,"class",sep="-")
      if ( is.null(res[[docname]]) ){
        res[[docname]] <- S4class.docs
        doc.names <- c(doc.names,docname)
      } else {
        stop(nn," appears as both S4 class and some other definition")
      }
    }
  }
  inherit.docs <- function(on){
    in.res <- res[[on]]
    if ( !is.null(parsed[[on]]) ){
      for ( parent in parsed[[on]]@parent ){
        if ( !is.na(parent) ){
          if ( is.null(in.res) ){
            in.res <- res[[parent]]
          } else if ( parent %in% names(res) ){
            parent.docs <- res[[parent]]
            for ( nn in names(parent.docs) ){
              if ( !nn %in% names(in.res) ){
                in.res[[nn]] <- parent.docs[[nn]]
              }
            }
          }
        }
      }
    }
    invisible(in.res)
  }
  all.done <- FALSE
  while ( !all.done ){
    res1 <- sapply(doc.names,inherit.docs,simplify=FALSE)
    all.done <- identical(res1,res)
    res <- res1
  }
  ## now strip out any generics (which have value NULL in res):
  res.not.null <- sapply(res,function(x){!is.null(x)})
  if ( 0 < length(res.not.null) && length(res.not.null) < length(res) ){
    res <- res[res.not.null]
  }
  res
### named list of lists, one for each object to document.
}

### List of parsers to use by default with package.skeleton.dx.
default.parsers <-
  c(extra.code.docs=extra.code.docs, ## TODO: cleanup!
    sapply(forfun.parsers,forfun),
    edit.package.file=function(desc,...){
      in.details <- setdiff(colnames(desc),"Description")
      details <- sprintf("%s: \\tab %s\\cr",in.details,desc[,in.details])
      L <-
        list(list(title=desc[,"Title"],
                  description=desc[,"Description"],
                  `tabular{ll}`=details))
      names(L) <- paste(desc[,"Package"],"-package",sep="")
      L
    },
    sapply(forall.parsers,forall)
    )

setClass("DocLink", # Link documentation among related functions
### The \code{.DocLink} class provides the basis for hooking together
### documentation of related classes/functions/objects. The aim is that
### documentation sections missing from the child are inherited from
### the parent class.
         representation(name="character", ##<< name of object
                        created="character", ##<< how created
                        parent="character", ##<< parent class or NA
                        code="character", ##<< actual source lines
                        description="character") ##<< preceding description block
         )

extract.file.parse <- function # File content analysis
### Using the base \code{parse} function, analyse the file to link
### preceding "prefix" comments to each active chunk. Those comments form
### the default description for that chunk. The analysis also looks for
### S4 class "setClass" calls and R.oo setConstructorS3 and setMethodS3
### calls in order to link the documentation of those properly.
(code
### Lines of R source code in a character vector - note that any
### nested \code{source} statements are \emph{ignored} when scanning
### for class definitions.
 ){
  res <- list()
  old.opt <- options(keep.source=TRUE)
  parsed <- try(parse(text=code))
  options(old.opt)
  if ( inherits(parsed,"try-error") ){
    stop("parse failed with error:\n",parsed)
  }
  chunks <- attr(parsed,"srcref")
  last.end <- 0
  for ( k in 1:length(parsed) ){
    start <- chunks[[k]][1]
    ##details<< If the definition chunk does not contain a
    ## description, any immediately preceding sequence consecutive
    ## "prefix" lines will be used instead.
    default.description <- NULL
    while ( start > last.end+1
           && grepl(prefix,code[start-1],perl=TRUE) ){
      start <- start-1
    }
    if ( start < chunks[[k]][1] ){
      default.description <- decomment(code[start:(chunks[[k]][1]-1)])
    } else {
      default.description <- NA_character_;
    }
    ##details<< Class and method definitions can take several forms,
    ## determined by expression type: \describe{
    ## \item{assignment (<-)}{Ordinary assignment of value/function;}
    ## \item{setClass}{Definition of S4 class;}
    ## \item{setConstructorS3}{Definition of S3 class using R.oo package;}
    ## \item{setMethodS3}{Definition of method for S3 class using R.oo package.}}
    ## Additionally, the value may be a name of a function defined elsewhere,
    ## in which case the documentation should be copied from that other definition.
    ## This is handled using the concept of documentation links.
    lang <- parsed[[k]]
    chars <- as.character(lang)
    expr.type <- chars[1]
    parent <- NA_character_

    if ( expr.type == "<-" || expr.type == "setConstructorS3" || expr.type == "setClass" ){
      object.name <- chars[2]
      ## If the function definition is not embedded within the call, then
      ## the parent is that function. Test whether the the third value
      ## looks like a name and add it to parents if so.
      if ( grepl("^[\\._\\w]+$",chars[3],perl=TRUE) ){
        parent <- chars[3]
      }
      res[[object.name]] <- new("DocLink",name=object.name,
                                created=expr.type,
                                parent=parent,
                                code=paste(chunks[[k]],sep=""),
                                description=default.description)
    } else if ( expr.type == "setMethodS3" || expr.type ==  "R.methodsS3::setMethodS3"){
      ##details<< The \code{setMethodS3} calls introduce additional
      ## complexity: they will define an additional S3 generic (which
      ## needs documentation to avoid warnings at package build time)
      ## unless one already exists. This also is handled by "linking"
      ## documentation. A previously unseen generic is linked to the
      ## first defining instances, subsequent definitions of that generic
      ## also link back to the first defining instance.
      generic.name <- chars[2]
      object.name <- paste(generic.name,chars[3],sep=".")
      if ( is.null(res[[generic.name]]) ){
        ## TDH 9 April 2012 Do NOT add \\link in generic.desc below,
        ## since it causes problems on R CMD check.
        ##* checking Rd cross-references ... WARNING
        ##Error in find.package(package, lib.loc) : 
        ##  there is no package called ‘MASS’
        ##Calls: <Anonymous> -> lapply -> FUN -> find.package

        generic.desc <-
          paste("Generic method behind \\code{",object.name,"}",sep="")
        res[[generic.name]] <- new("DocLink",
                                   name=generic.name,
                                   created=expr.type,
                                   parent=object.name,
                                   code=NA_character_,
                                   description=generic.desc)
      } else {
        parent <- res[[generic.name]]@parent
      }
      ## If the function definition is not embedded within the call, then
      ## the parent is that function. Test whether the the fourth value
      ## looks like a name and add it to parents if so.
      if ( grepl("^[\\._\\w]+$",chars[4],perl=TRUE) ){
        parent <- c(chars[4],parent)
      }
      res[[object.name]] <- new("DocLink",name=object.name,
                                created=expr.type,
                                parent=parent,
                                code=paste(chunks[[k]],sep=""),
                                description=default.description)
    } else {
      ## Not sure what to do with these yet. Need to deal with setMethod, setAs etc.
    }
  }
  invisible(res)
### Returns an invisible list of .DocLink objects.
}

extract.docs.setClass <- function # S4 class inline documentation
### Using the same conventions as for functions, definitions of S4 classes
### in the form \code{setClass("classname",\dots)} are also located and
### scanned for inline comments.
(doc.link
### DocLink object as created by \code{extract.file.parse}.
### Note that \code{source} statements are \emph{ignored} when scanning for
### class definitions.
 ){
  chunk.source <- doc.link@code
  ##details<<
  ## Extraction of S4 class documentation is currently limited to expressions
  ## within the source code which have first line starting with
  ## \code{setClass("classname"}. These are located from the source file
  ## (allowing also for white space around the \code{setClass} and \code{(}).
  ## Note that \code{"classname"} must be a quoted character string;
  ## expressions returning such a string are not matched.
  class.name <- doc.link@name

  ##details<< For class definitions, the slots (elements of the
  ## \code{representation} list) fill the role of function
  ## arguments, so may be documented by \code{##<<} comments on
  ## the same line or \code{### } comments at the beginning of the
  ## following line.
  f.n <- paste(class.name,"class",sep="-")
  docs <- extract.xxx.chunks(chunk.source,f.n)
  ## also apply source parsing functions that I separated out into
  ## separate functions
  docs <- combine(docs,lonely$prefixed.lines(chunk.source))
  docs$title <- lonely$title.from.firstline(chunk.source)
  ##details<<
  ## If there is no explicit title on the first line of setClass, then
  ## one is made up from the class name.
  if ( 0 == length(docs$title) ){
    docs$title <- list(title=paste(class.name,"S4 class"))
  }
  ##details<<
  ## The class definition skeleton includes an \code{Objects from the Class}
  ## section, to which any \code{##details<<} documentation chunks are
  ## written. It is given a vanilla content if there are no specific
  ## \code{##details<<} documentation chunks.
  if ( is.null(docs[["details"]]) ){
    docs[["details"]] <-
      paste("Objects can be created by calls of the form \\code{new(",
            class.name," ...)}",sep="")
  }
  docs[["section{Objects from the Class}"]] <- docs[["details"]]
  ## seealso has a skeleton line not marked by ~ .. ~, so have to suppress
  if ( is.null(docs[["seealso"]]) ){
    docs[["seealso"]] <- ""
  }
  if ( is.null(docs[["alias"]]) ){
    docs[["alias"]] <- class.name
  }
  if ( is.null(docs[["description"]]) ){
    docs[["description"]] <- doc.link@description
  }
  invisible(docs)
}

apply.parsers <- function
### Parse code to r objs, then run all the parsers and return the
### documentation list.
(code,
### Character vector of code lines.
 parsers=default.parsers,
### List of Parser Functions.
 verbose=FALSE,
### Echo names of Parser Functions?
 ...
### Additional arguments to pass to Parser Functions.
 ){
  e <- new.env()
  ## KMP 2011-03-09 fix problem with DocLink when inlinedocs ran on itself
  ## Error in assignClassDef(Class, classDef, where) :
  ##   Class "DocLink" has a locked definition in package "inlinedocs"
  ## Traced to "where" argument in setClassDef which defaults to topenv()
  ## which in turn is inlinedocs when processing inlinedocs package, hence
  ## the clash. The following works (under R 2.12.2), so that the topenv()
  ## now finds e before finding the inlinedocs environment.
  old <- options(keep.source=TRUE,topLevelEnvironment=e)
  on.exit(options(old))
  exprs <- parse(text=code)
  ## TDH 2011-04-07 set this so that no warnings about creating a fake
  ## package when we try to process S4 classes defined in code
  e$.packageName <- "inlinedocs.processor"
  for (i in exprs){
      eval(i, e)
  }
  objs <- sapply(ls(e),get,e,simplify=FALSE)

  docs <- list()

  ## apply parsers in sequence to code and objs
  if(verbose)cat("Applying parsers:\n")
  for(i in seq_along(parsers)){
    N <- names(parsers[i])
    #mm if(verbose){
      if(is.character(N) && N!=""){
        cat(" this is parser:",N,"\n",sep="")
      }else cat('.\n')
    #mm }
    p <- parsers[[i]]
    ## This is the argument list that each parser receives:
    L <- p(code=code,objs=objs,docs=docs,env=e,...)
    # print("mm point1")
    #save(docs,L,file="/home/mm/SoilR/scripts/docs_L.RData")
    #print(paste(L,"\n"))
    #if(N=="exclude")browser()
    docs <- combine(docs,L) #mm
  }
  ## post-process to collapse all character vectors
  for(i in seq_along(docs)){
    for(j in seq_along(docs[[i]])){
      if(names(docs[[i]])[j]!=".s3method")
      docs[[i]][[j]] <- paste(docs[[i]][[j]],collapse="\n")
    }
 }
  if(verbose)cat("\n")
  return(docs)
### A list of extracted documentation from code.
}

### Names of Parser Functions that operate on the desc arg.
descfile.names <- c("author.from.description","edit.package.file")

### Names of Parser Functions that do NOT use the desc arg.
non.descfile.names <-
  names(default.parsers)[!names(default.parsers)%in%descfile.names]

### Parsers that operate only on R code, independently of the
### description file.
nondesc.parsers <- default.parsers[non.descfile.names]

extract.docs.file <- structure(function
### Apply all parsers relevant to extract info from just 1 code file.
(f,
### File name of R code to read and parse.
 parsers=NULL,
### Parser Functions to use to parse the code and extract
### documentation.
 ...
### Other arguments to pass to Parser Functions.
 ){
  if(is.null(parsers))parsers <- nondesc.parsers
  apply.parsers(readLines(f),parsers,verbose=FALSE,...)
},ex=function(){
  f <- system.file("silly","R","silly.R",package="inlinedocs")
  extract.docs.file(f)
})

