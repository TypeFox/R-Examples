package.skeleton.dx <- structure(function # Package skeleton deluxe
### Generates Rd files for a package based on R code and DESCRIPTION
### metadata. After inspecting the specified R code files to find
### inline documentation, it calls the standard package.skeleton
### function, which creates bare Rd files. The inline documentation is
### added to these Rd files and then these files are copied to
### pkgdir/man, possibly overwriting the previous files there.
(pkgdir="..",
### Package directory where the DESCRIPTION file lives. Your code
### should be in pkgdir/R. We will setwd to pkgdir/R for the duration
### of the function, then switch back to where you were previously.
 parsers=NULL,
### List of Parser functions, which will be applied in sequence to
### extract documentation from your code. Default NULL means to first
### search for a definition in the variable "parsers" in
### pkgdir/R/.inlinedocs.R, if that file exists. If not, we use the
### list defined in options("inlinedocs.parsers"), if that is
### defined. If not, we use the package global default in the variable
### default.parsers.
 namespace = FALSE,
### A logical indicating whether a NAMESPACE file should be generated
### for this package. If \code{TRUE}, all objects whose name starts
### with a letter, plus all S4 methods and classes are exported.
 excludePattern=FALSE,
### A regular expression matching the files that are not to be
### processed e.g. because inlinedocs can not handle them yet (like
### generic function definitions)
 ...
### Parameters to pass to Parser Functions.
 ){
  ## This causes a warning on R CMD check TDH 28 Jan 2013.
  ##alias< < inlinedocs
	 
  chdir <- file.path(pkgdir,"R")
  if(!file.exists(chdir))stop("need pkgdir/R, tried ",chdir)
  old.wd <- setwd(chdir)
  on.exit(setwd(old.wd))
  # PhG: R allows for specific code to be in /unix, or /windows subdirectories
  # but apparently, inlinedocs does not support this. I think it is fair to
  # stop here with an explicit error message if at least one of /unix or
  # /windows subdirectory is found!
  # file_test(-d, ...) does the job, but I don't want to add a dependency on
  # package 'utils", where it lives. So, I prefer using file.info()
  if (isTRUE(file.info("unix")$isdir) || isTRUE(file.info("windows")$isdir))
    stop("Platform-specific code in ./R/unix, or ./R/windows is not supported")

  ## Default values and required fields in DESCRIPTION file.
  description.defaults <-
    c("Package"="",
      "Maintainer"=Sys.getenv("USER"),
      "Author"=Sys.getenv("USER"),
      "Version"="1.0",
      "License"="GPL-3",
      "Title"="a package",
      "Description"="a package that does\n many things.")

  ## Necessary fields in DESCRIPTION, otherwise error.
  fields <- names(description.defaults)

  ## Default DESCRIPTION, written if it doesn't exist.
  empty.description <-
    matrix(description.defaults,ncol=length(fields),dimnames=list(NULL,fields))


  
  
  ## if no DESCRIPTION, make one and exit.
  descfile <- file.path("..","DESCRIPTION")
  if(!file.exists(descfile)){
    write.dcf(empty.description,descfile)
    stop("Need ",descfile,"; please fill that in and try again")
  }

  ## Read description and check for errors
  desc <- read.dcf(descfile)
  ## TDH 3 Sept 2013 need to support Authors@R for CRAN.
  if("Authors@R" %in% colnames(desc)){
    author <- paste(eval(parse(text=desc[,"Authors@R"])), collapse=", ")
    desc <- cbind(desc,
                  Author=author,
                  Maintainer=author)
  }
  if(any(f <- !sapply(fields,is.element,colnames(desc))))
    stop("Need ", paste(names(f)[f], collapse = ", "), " in ", descfile)
    #PhG: corrected from stop("Need ",names(f)[f]," in ",descfile)
  if(any(f <- sapply(fields,function(f)desc[,f]=="")))
    stop("Need a value for ", paste(names(f)[f], collapse = ", "),
         " in ", descfile)
    #PhG: corrected from stop("Need a value for ",names(f)[f]," in ",descfile)

  ## Load necessary packages before loading pkg code
  if("Depends" %in% colnames(desc)){
    required <- strsplit(desc[,"Depends"],split=",")[[1]]
    # PhG: for packages with NAMESPACE, dependencies are also listes in the
    # Imports field!
  } else required <- character(0)
  # PhG: packages listed in Imports field are not supposed to be attached to the
  # search path when a package with NAMESPACE is loaded, only the namespace is
  # loaded. However, the code here is not loaded as it should be, and we need to
  # load also these Import(ed) packages to get correct results in the present
  # case (to be checked with most complex cases!)
  if ("Imports" %in% colnames(desc))
    required <- c(required, strsplit(desc[, "Imports"], split = ",")[[1]])
  ##print(c("required=",required))
  ## This may need some refining, basically I just tried to take the
  ## first word from each vector element, stripping of whitespace in
  ## front and anything after:
  #pkgnames <- gsub("\\W*(\\w+)\\b.*","\\1",required)
  # PhG: the previous line is wrong: it does not work with package names
  # like R.oo... Extract from Writing R Extensions manual:
  # "The `Package' and `Version' fields give the name and the version of the
  # package, respectively. The name should consist of letters, numbers, and the
  # dot character and start with a letter." 
  # Consequently, I propose:
  pkgnames <- gsub("\\W*([a-zA-Z][a-zA-Z0-9.]*)\\b.*", "\\1", required)
  # PhG: We need to eliminate 'R' from the list!
  pkgnames <- pkgnames[pkgnames != "R"]
  # PhG: if we create a namespace, we need to keep this list for further use
  if (isTRUE(namespace)) allpkgs <- pkgnames
  # PhG: We eliminate also from the list the packages that are already loaded
  if(length(pkgnames))
    pkgnames <- pkgnames[!sprintf("package:%s",pkgnames) %in% search()]
  # PhG: according to Writing R Extensions manual, a package name can occur
  # several times in Depends
  pkgnames <- unique(pkgnames)
  if (length(pkgnames)) {
    # PhG: A civilized function returns the system in the same state it was
    # before => detach loaded packages at the end!
    on.exit(suppressWarnings({
      try(for (pkg in pkgnames)detach(paste("package", pkg, sep = ":"),
        unload = TRUE, character.only = TRUE), silent = TRUE)
    }), add = TRUE)
    # PhG: Shouldn't we need to check that packages are loaded and shouldn't
    # we exit with an explicit error message if not? Note: we don't use version
    # information here. That means we may well load wrong version of the
    # packages... and that is NOT detected as an error!
    pkgmissing <- character(0)
    for (pkg in pkgnames) {
      if(!require(pkg, character.only = TRUE)){
        pkgmissing <- c(pkgmissing, pkg)
      }
    }
    if (length(pkgmissing))
        stop("Need missing package(s): ", paste(pkgmissing, collapse = ", "))
  }

  ## for the parser list, first try reading package-specific
  ## configuration file
  if(is.null(parsers))parsers <- tryCatch({
    cfg <- new.env()
    sys.source(cfile <- ".inlinedocs.R",cfg)
    L <- cfg$parsers
    if(!is.null(L))cat("Using parsers in ",cfile,"\n",sep="")
    L
  },error=function(e)NULL)
  ## then try the global options()
  opt <- "inlinedocs.parsers"
  if(is.null(parsers)&&!is.null(parsers <- getOption(opt))){
    cat("Using parsers in option ",opt,"\n",sep="")
  }
  ## if nothing configured, just use the pkg default
  if(is.null(parsers))parsers <- default.parsers
  
  ## concatenate code files and parse them
  # PhG: in Writing R Extensions manuals, source code in /R subdirectory can
  # have .R, .S, .q, .r, or .s extension. However, it makes sense to restrict
  # this to .R only for inlinedocs, but a clear indication is required in the
  # man page!
  code_files <- if(!"Collate"%in%colnames(desc))Sys.glob("*.R")
  else strsplit(gsub("\\s+"," ",desc[,"Collate"]),split=" ")[[1]]
  code_files =grep(excludePattern,code_files,invert=TRUE,value=TRUE)
  ## TDH 28 Jan 2013, warn users such as Pierre Neuvial if they have
  ## comments on the last line of one input file. Sometimes comments
  ## on the last line can appear to be the first line of comments of
  ## the next code file.
  lines.list <- lapply(code_files,readLines)
  for(i in seq_along(lines.list)){
    lvec <- lines.list[[i]]
    fn <- code_files[i]
    last <- lvec[length(lvec)]
    if(grepl("#",last)){
      warning("comment on last line of ",fn,
              ", unexpected docs may be extracted")
    }
  }
  #print(excludePattern)
  ## Make package skeleton and edit Rd files (eventually just don't
  ## use package.skeleton at all?)
  name <- desc[,"Package"]
  unlink(name,recursive=TRUE)
  package.skeleton(name,code_files=code_files)

#  # PhG: one must consider a potential Encoding field in DESCRIPTION file!
  # which is used also for .R files according to Writing R Extensions
  if ("Encoding" %in% colnames(desc)) {
    oEnc <- options(encoding = desc[1, "Encoding"])$encoding
    on.exit(options(encoding = oEnc), add = TRUE)
  }
  code <- do.call(c,lapply(code_files,readLines))
#  #print(code)
  docs<- apply.parsers(code,parsers,verbose=TRUE,desc=desc)

  cat("Modifying files automatically generated by package.skeleton:\n")
  ## documentation of generics may be duplicated among source files.
  dup.names <- duplicated(names(docs))
  if ( any(dup.names) ){
    warning("duplicated file names in docs: ",paste(names(docs)[dup.names]))
  }
  for(N in unique(names(docs))) {
	modify.Rd.file(N,name,docs)
        if(grepl("class",N)){
		removeAliasesfrom.Rd.file(N,name,code)
	}
  }
  file.copy(file.path(name,'man'),"..",recursive=TRUE)
  # PhG: copy NAMESPACE file back
  if (isTRUE(namespace)) {
    # PhG: package.skeleton() does not add import() statement, but here the
    # philosophy is to get a fully compilable package, which is not the
    # case at this stage with a NAMESPACE. So, we add all packages listed
    # in Depends and Imports fields of the DESCRIPTION file in an import()
    # statement in the NAMESPACE
    nmspFile <- file.path("..", "NAMESPACE")
    cat("import(", paste(allpkgs, collapse = ", "), ")\n\n", sep = "",
        file = nmspFile) 
    # PhG: append the content of the NAMESPACE file generated by
    # package.skeleton()
    file.append(nmspFile, file.path(name,'NAMESPACE'))
    # PhG: we also have to export S3 methods explictly in the NAMESPACE
    cat("\n", file = nmspFile, append = TRUE)
    for (N in unique(names(docs))) {
        d <- docs[[N]]
        if (!is.null(d$.s3method))
            cat('S3method("', d$.s3method[1], '", "', d$.s3method[2], '")\n',
                sep = "", file = nmspFile, append = TRUE)
    }
  }
  
  unlink(name,recursive=TRUE)
},ex=function(){
  owd <- setwd(tempdir())
  
  ## get the path to the silly example package that is provided with
  ## package inlinedocs
  testPackagePath <- file.path(system.file(package="inlinedocs"),"silly")
  ## copy example project to the current unlocked workspace that can
  ## be modified
  file.copy(testPackagePath,".",recursive=TRUE)
  
  ## generate documentation .Rd files for this package
  package.skeleton.dx("silly")
  
  ## check the package to see if generated documentation passes
  ## without WARNINGs.
  if(interactive()){
    cmd <- sprintf("%s CMD check --as-cran silly",file.path(R.home("bin"), "R"))
    print(cmd)
    checkLines <- system(cmd,intern=TRUE)
    warnLines <- grep("WARNING",checkLines,value=TRUE)
    if(length(warnLines)>0){
      writeLines(checkLines)
      cat("\n\nLines with WARNING:\n")
      print(warnLines)
      ## disable due to bug in R CMD check:
      ## https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14875
      ##stop("WARNING encountered in package check!")
    }
  }  
  ## cleanup: remove the test package from current workspace again
  unlink("silly",recursive=TRUE)
  setwd(owd)
})


replace.one <- function
### Do find and replace for one element of an inner documentation list
### on 1 Rd file.
(torep,
### tag to find.
 REP,
### contents of tag to put inside.
 txt
### text in which to search.
 ){
  ##if(grepl("Using the same conventions",REP))browser()
  escape.backslashes <- function(x)gsub("\\\\","\\\\\\\\",x)
  cat(" ",torep,sep="")
  FIND1 <- escape.backslashes(torep)
  FIND <- gsub("([{}])","\\\\\\1",FIND1)
  FIND <- paste(FIND,"[{][^}]*[}]",sep="")
  REP.esc <- escape.backslashes(REP)
  ## need to escape backslashes for faithful copying of the comments
  ## to the Rd file:
  REP <- paste(FIND1,"{",REP.esc,"}",sep="")
  ## escape percent signs in R code:
  REP <- gsub("%","\\\\\\\\%",REP)
  ## alias (in particular) need to change only the first one generated
  ## (generic methods in classes add to standard skeleton alias set)
  if ( torep %in% c("alias") ){
    txt <- sub(FIND,REP,txt)
  } else {
    txt <- gsub(FIND,REP,txt)
  }
  classrep <- sub("item{(.*)}","item{\\\\code{\\1}:}",torep,perl=TRUE)
  if ( classrep != torep ){
    ## in xxx-class files, slots are documented with:
    ## \item{\code{name}:}{Object of class \code{"function"} ~~ }
    ## which requires slightly different processing
    FIND1 <- escape.backslashes(classrep)
    FIND <-
      paste(gsub("([{}])","\\\\\\1",FIND1),
            "\\{Object of class \\\\code\\{\\\"(\\S+)\\\"\\}[^}]*[}]",sep="")
    ## need to escape backslashes for faithful copying of the comments
    ## to the Rd file and also put the class type in parentheses.
    REP <- paste(FIND1,"{(\\\\code{\\1}) ",REP.esc,"}",sep="")
    ## escape percent signs in R code:
    REP <- gsub("%","\\\\\\\\%",REP)
    txt <- gsub(FIND,REP,txt)
  }
  txt
}

removeAliasesfrom.Rd.file <- function
### remove aliases to methodnames from the Rd file of a class
### automatically-generated by package.skeleton.
(N,
### Name of function/file to which we will add documentation.
 pkg,
### Package name.
code
### The code of the package 
)
{
  # <mm:package.skeleton adds some duplicated aliases to the .*-class.Rd files 
  # to get rid of the warnings from R CMD check 
  # we have to delete them 
  # the duplicated aliases have 2 sources

  Nme <- fixPackageFileNames(N)
  fb <- paste(Nme,".Rd",sep="")
  f <- file.path(pkg,'man',fb)

  ## If we do.not.generate this file, it does not exist so we need to
  ## do nothing to avoid errors.
  if(!file.exists(f))return()
  
  dlines <- readLines(f)
  
  name=gsub("-class","",N)
  
  # these we will later comment out
  # now we look at all the aliases produces by package.skeleton
  aliasLine2name=function(line){return(gsub("^\\\\alias\\{(.*)\\}","\\1",line))}
  aliasInd <- grep("^\\\\alias.*",dlines)
  # aliasnames=as.character(lapply(dlines[aliasInd],aliasLine2name))
  # duplicates= intersect(methodnames,aliasnames)
  
  if(length(aliasInd)){
      # first get rid of the ${classname}-mehthod stuff
      p1=paste(",",name,"-method",sep="")
  
      # next get rid of the -mehthod stuff
      patterns=c(p1,",ANY-method",",ANY",",character")
      for (pattern in patterns){
          dlines[aliasInd] <- gsub(pattern,"",dlines[aliasInd]) 
      }
  }
  # next: 
  # the names of the methods implemented by the class
  # which occur in all ".*-class.Rd" files 
  # as aliases 
  # to find them we ask for those methods but 
  # therefore have to readthe the code to be documented.
  # As the apply.parsers function we do this in a separate 
  # environment
  # mm:
  # This could probably be factored out since it duplicates
  # some functionality of apply.parsers (evaluating the sources again)
  # but I am not sure how big the changes involved woud be.
  
  e <- new.env()
  old <- options(keep.source=TRUE,topLevelEnvironment=e)
  on.exit(options(old))
  exprs <- parse(text=code)
  ## set this so that no warnings about creating a fake
  ## package when we try to process S4 classes defined in code
  e$.packageName <- "inlinedocs.processor"
  for (i in exprs){
      eval(i, e)
  }
  g=file(,open="w+")#anonymous file
  showMethods(classes=name,printTo=g,where=e)
  lines=readLines(g)
  close(g)
  ind=grep("Function",lines)
  flines=lines[ind]
  line2name=function(line){return(strsplit(line," ")[[1]][2])}
  methodnames=as.character(lapply(flines,line2name))
  markDupAliases=function(line){
  	ret=line
  	if(is.element(aliasLine2name(line),methodnames)){
  		ret=gsub("(^\\\\alias)","%% \\1",line)
  	}
  	return(ret)
  } 
  for (j in aliasInd){
  	dlines[[j]]=markDupAliases(dlines[[j]])
  }
  #print(fb)
  #print(d)
  #print(dlines)
  txt <- paste(dlines,collapse="\n")
  fc=file(f,open="w+")#anonymous file
  writeLines(txt,fc)
  close(fc)
}	

modify.Rd.file <- function
### Add inline documentation from comments to an Rd file
### automatically-generated by package.skeleton.
(N,
### Name of function/file to which we will add documentation.
 pkg,
### Package name.
 docs
### Named list of documentation in extracted comments.
 ){
  # PhG: for functions like 'obj<-', package.skeleton creates files like 'obj_-'
  # => rework names the same way, i.e., using the same function from utils package
  Nme <- fixPackageFileNames(N)
  fb <- paste(Nme,".Rd",sep="")
  ## For some functions, such as `[[.object`, package.skeleton (as used
  ## within this package but not when used standalone) seems to generate
  ## with a preceding z ("z[[.object.Rd"), so the z form is tested for and
  ## used if it exists and the first does not.
  zfb <- paste("z",Nme,".Rd",sep="")
  f <- file.path(pkg,'man',fb)
  if ( (!file.exists(f)) && file.exists(file.path(pkg,'man',zfb)) ){
    fb <- zfb
    f <- file.path(pkg,'man',zfb)
  }
  ## If there are no significant docs in the comments then the object
  ## should still be documented, by writing the file by hand in the
  ## man directory. This will write a blank Rd file if none exists, so
  ## it's easy to get started.
  if((length(docs[[N]])<3) &&
     file.exists(file.path("..","man",fb))){
    # print(paste("mm object with no documentation available N=",N))
    unlink(f)
    return()
  }
  cat(N,":",sep="")
  d <- docs[[N]]
  ## for some functions no documentatian file is created by package.skeleton
  ## for instance generic functions that are already defined in other packages
  ## like print or plot so there is still the possibillity that 
  ## f is missing altogether
  if (!file.exists(f)) {
	print("mm missing file")
	print(f)
	return()
	}
  dlines <- readLines(f)

  ## cut out alias line if we are in the package file and there is a
  ## matching function
  if(length(grep("-package$",N)) && "alias" %in% names(d) )
    dlines <- dlines[-grep(paste("alias[{]",N,sep=""),dlines)-1]
  else if ( "alias" %in% names(d) ){
    ## allowing alias changes have to make sure that original alias remains
    ## note that the contents of this go inside \alias{}, so the separator
    ## has to terminate one and start the next
    d[["alias"]] <- paste(paste(N,"}\n\\alias{",sep=""),
                            d[["alias"]],sep="")
  }
  

  # PhG: in the special case of custom operators like %....%, we must protect
  # these strings in name, alias and usage (at least)! Otherwise, bad things
  # happen with these strings: (1) usage entry is cut out, because confused
  # with comments, and % are escaped in name and alias!
  if (grepl("^%.+%$", N)) {
    Nmask <- gsub("%", "---percent---", N)
    # Replace any occurence of N by Nmask
    dlines <- gsub(N, Nmask, dlines, fixed = TRUE)
  } else Nmask <- NULL

  ## cut out all comments {} interferes with regex matching
  comments <- grep("^[%~]",dlines)
  ## gotcha! if no comment lines, then -(nothing) kills everything
  if ( 0 < length(comments) ) dlines <- dlines[-comments]
  ## and class skeletons have a different way of using ~
  dlines <- gsub("\\s*~.*~\\s*","",dlines,perl=TRUE)
  ## and the "Objects can be created..." boilerplate also breaks perl REs
  dlines <- gsub("Objects can be created by calls.*\\)\\}","",dlines,perl=TRUE)
  ## ditto the "or \code{\linkS4class{CLASSNAME}} for links to other classes"
  dlines <- gsub("or \\\\code\\{\\\\linkS4class.*classes","",dlines,perl=TRUE)

  ## cut out a couple of sections that cause warnings
  o <- grep("Optionally",dlines)
  if(length(o))dlines <- dlines[-(o:(o+1))]
  ## delete examples til the end of the file (also includes keywords)
  dlines <- dlines[1:(tail(grep("examples[{]$",dlines),1)-1)]
  ## add back a minimal examples section to find and replace
  dlines <- c(dlines,"\\examples{}\n")
  ## and replace keyword section if keywords are present.
  if ( "keyword" %in% names(d) ){
    dlines <- c(dlines,"\\keyword{}\n")
  }

  ## erase curly braces in format section, which appear sporadically
  ## and can cause errors in R CMD check.
  fstart <- grep("^\\\\format[{]$",dlines)+1
  if(length(fstart)){
    closing <- grep("^[}]$",dlines)
    fend <- closing[closing>fstart][1]-1
    dlines[fstart:fend] <- gsub("[{}]","",dlines[fstart:fend])
  }

  ## sometimes (s4 classes) title is has \code{} blocks inside, which
  ## causes problems with our find-replace regex inside replace.one,
  ## so lets just put a simple title that works.
  i <- grep("^\\\\title",dlines)
  if(length(i)){
    dlines[i] <- gsub("\\\\code[{][^}]*[}]","",dlines[i])
  }
  name=N
  txt <- paste(dlines,collapse="\n")
  
  ## Fix usage
  m <- regexpr("usage[{][^}]*[}]",txt)
  Mend <- m+attr(m,"match.length")
  utxt <- substr(txt,m+6,Mend-2)

  ## fix \method version if .s3method
  if ( is.null(d$.s3method) ) {
	  # PhG: in case we have fun<-(x, ..., value), we must rewrite it
	  # as fun(x, ...) <- value
	  if (grepl("<-$", N)) {
		  utxt <- sub("<-[(](.+), ([^,)]+)[)]",
				  "(\\1) <- \\2", utxt)
	  }
	  # PhG: this is for special functions %...% which should write x %...% y
	  if (grepl("^%.*%$", N)) {
		  utxt <- sub("(%.*%)[(]([^,]+), ([^)]+)[)]",
				  "\\2 \\1 \\3", utxt) 
	  }
  }
  
  ## multiple lines for the PDF!
  # tw: parse fails on accessor functions such as "myF<-" <- function(data,x) 
  # see testfile accessorFunctions.R
  # workaround with tryCatch
  parsed <- utxt
  tryCatch({
	  parsed <- parse(text=utxt)
  }, error = function(e) warning(e) )
  if(length(parsed)){
	  utxt <- sprintf("usage{%s}\n",paste(format(parsed[[1]]),collapse="\n"))
  }
  if(length(grep("usage[{]data",utxt))){
    utxt <- gsub("data[(]([^)]*)[)]","\\1",utxt)
  }
  ## fix \method version if .s3method

  if ( !is.null(d$.s3method) ){
    pat <- paste(d$.s3method,collapse=".")
    rep <- paste("\\method{xx",d$.s3method[1],"}{",d$.s3method[2],"}",sep="")
    utxt <- gsub(pat,rep,utxt,fixed=TRUE)
    
    # PhG: there is the special case of generic<-.obj(x, ..., value) to rewrite
    # \method{generic}{obj}(x, ...) <- value
    if (grepl("<-$", d$.s3method[1])) {
        # 1) replace {generic<-} by {generic}
        utxt <- sub("<-[}]", "}", utxt)
        # 2) replace ..., value) by ...) <- value
        utxt <- sub(", *([^),]+)[)]", ") <- \\1", utxt)
    }
  } else {
	#tw: moved before parse
  }
  ## add another backslash due to bug in package.skeleton
  ## but only if not before % character due to another bug if % in usage
  ## arguments - see above
  txt <- paste(substr(txt,1,m-1),
               gsub("\\\\([^%])","\\\\\\\\\\1",utxt),
               substr(txt,Mend+1,nchar(txt)),
               sep="")


  
  ## At least in my code, any remaining % symbols are in \usage sections
  ## as function arguments. These promptly break Rd check because you end
  ## up with unterminated strings. Just in case, the following regexp only
  ## modifies those % symbols which follow something other than %.
  ## (a more complicated version would attempt to do so only within strings.)
  txt <- gsub("([^%])%","\\1\\\\%",txt,perl=TRUE)
  
  # PhG: now restore masked function name, if any (case of %....% operators)
  if (!is.null(Nmask))
    txt <- gsub(Nmask, N, txt, fixed = TRUE)

  ## Find and replace based on data in d
  for(torep in names(d)){
    if ( !grepl("^[.]",torep) ){## .flags should not be used for find-replace
      txt <- replace.one(torep,d[[torep]],txt)
    }
  }

  ## delete empty sections to suppress warnings in R CMD check
  txt <- gsub("\\\\[a-z]+[{]\\s*[}]","",txt)
  if ( !is.null(d$.s3method) ){
    ## and now remove the xx inserted above to prevent \method{[[}{...} falling
    ## foul of the above replacement!
    txt <- gsub("\\\\method{xx","\\method{",txt,fixed=TRUE)
  }
  ## This doesn't work if there are quotes in the default values:
  ## gsub(",",paste("\n",paste(rep(" ",l=nchar(N)-1),collapse="")),utxt)

  ## convert to dos line endings to avoid problems with svn
  txt <- gsub("(?<!\r)\n","\r\n",txt,perl=TRUE)
  cat(txt,file=f)
  cat("\n")
}
