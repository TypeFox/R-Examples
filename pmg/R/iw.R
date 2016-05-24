##file import Wizard
##Uses BasicWidgets
##call pmg.specifyFileForImport to start it off

.l = list()
.l[[gettext("text files")]] = c("csv","txt","fwf")
.l[[gettext("ARFF files")]] = c("arff")
.l[[gettext("DBF files")]] = c("dbf")
.l[[gettext("Stata Binary files")]] = c("dta")
.l[[gettext("EPI info files")]] = c("epi")
.l[[gettext("Minitab Portable files")]] = c("mtp")
.l[[gettext("Octave text data files")]] = c("octave")
.l[[gettext("SPSS files")]] = c("sav")
.l[[gettext("SAS XPORT files")]] = c("xport")
.l[[gettext("Systat files")]] = c("sys","syd")
.l[[gettext("Excel files")]] = c("xls")
.l[[gettext("DIF files")]] = c("DIF","dif")
.l[[gettext("Open office files")]] = c("odt")
.l[[gettext("gnumeric files")]] = c("gnumeric")
.fileExtensions =  .l

## .fileExtensions =  list(
##   "text files" = c("csv","txt","fwf"),
##   "ARFF files" = c("arff"),
##   "DBF files" = c("dbf"),
##   "Stata Binary files" = c("dta"),
##   "EPI info files" = c("epi"),
##   "Minitab Portable files" = c("mtp"),
##   "Octave text data files" = c("octave"),
##   "SPSS files" = c("sav"),
##   "SAS XPORT files" = c("xport"),
##   "Systat files" = c("sys","syd"),
##   "Excel files" = c("xls"),
##   "DIF files" = c("DIF","dif"),
##   "Open office files" = c("odt"),
##   "gnumeric files" = c("gnumeric")
##   )
## strip last character
pop = function(x) x[-length(x)]
popchar = function(str) paste(pop(unlist(strsplit(str,""))),collapse="")

selectFile = function(initialFile = NULL) {

  filterList = lapply(.fileExtensions, function(i) list(patterns = paste("*.",i,sep="")))
  filterList$"All files" = list(patterns=c("*"))
  gfile(text = "Select a file for import",
        initialfilename = initialFile,
        filter = filterList
        )
}


## specify with a URL or a filebrowse
pmg.specifyFileForImport = function(...) {

  filterList = lapply(.fileExtensions, function(i) list(patterns = paste("*.",i,sep="")))
  filterList$"All files" = list(patterns=c("*"))

  GUI = BasicGUI$new(message=gettext("Select a file to import"))
  GUI$filterList = filterList
  GUI$useDefaultText = gettext("<use file extension to determine>")
  GUI$fileSelectDefaultText = gettext("Specify a file or url...")
  GUI$makeBody = function(., container) {
    g = ggroup(horizontal=FALSE, container=container)
    tbl = glayout(container=g)
    tbl[1,1] <- "local file"
    tbl[1,2] <- (.$filebrowse = gfilebrowse(text=.$fileSelectDefaultText,
                   action=invisible,
                   container=tbl, filter=.$filterList, quote=FALSE))
    tbl[2,1] <- (l <- glabel(gettext("or"),container=tbl))
    font(l) <- c(style="italic")
    tbl[2,2] <- gseparator(container=tbl)
    tbl[3,1] <- "url"
    tbl[3,2] <- (.$url = gedit("", container=tbl))

    tbl[4,1:2] <- gseparator(container=tbl)
    tbl[5,1] = gettext("File type is")
    tbl[5,2] <- (.$filetype = gdroplist(c(
      .$useDefaultText,
      sapply(names(filterList),popchar)
      ), container=tbl))

    visible(tbl) <- TRUE
  }
  GUI$clearButtonHandler = NULL
  GUI$okButtonHandler = function(.,h,...) {
    ## what to do? need *local* filename and type
    ## if url, but no file, then we download file name it, go
    ## if file then go to next

    .$theFile = svalue(.$filebrowse)
    theURL = svalue(.$url)
    .$ext = NULL ## the extension, figure out


    if(.$theFile == .$fileSelectDefaultText || !file.exists(.$theFile)) {
      ## try to get the URL
      .$theFile= tempfile()
      out = try(download.file(theURL, destfile = .$theFile))
      if(inherits(out,"try-error")) {
        sprintf("Error downloading file: %s\n",out)
        return(TRUE)
      }
      ## we saved to out
      ## guess extension from $url
      tmp = unlist(strsplit(basename(theURL), split="\\."))
      .$ext = tmp[length(tmp)]
    }
    ##  file is now theFile
    ## get extension type from droplist

    fileType = svalue(.$filetype)

    if(fileType != .$useDefaultText) {
      ## use filterList to get
      fileType = paste(fileType,"s", sep="", collapse="") ## append s back
      .$ext = .fileExtensions[[fileType]][1]
      sprintf("Set extension to %s \n",.$ext)
    } else if(is.null(.$ext)) {
      tmp = unlist(strsplit(basename(.$theFile), split="\\."))
      .$ext = tmp[length(tmp)]
    } 
    ## now we have .$theFile and .$ext move on
    dispose(.$window)

    importFile(.$theFile, .$ext)
  }

  ## now draw GUI
  GUI$show()
}


importFile = function(filename, ext=NULL) {

  if(missing(filename))
    filename = selectFile()



  GUI = BasicGUI$new(message=paste("import", filename,collapse=" "))
  GUI$filename = filename
  GUI$ext = ext
  GUI$AssignToText = gettext("Assign to:")
  GUI$clearButtonHandler = NULL
  GUI$okButtonHandler = function(.,h,...) {
    ## the functions below define FUN, args, and varName
    out = try(do.call(.$FUN,lapply(args,svalue)), silent=TRUE)
    if(inherits(out,"try-error")) {
      sprintf("Error: %s \n",out)
    } else {
      varName = make.names(svalue(.$varName))
      ## can't have empty names due to make.names

      ## check if there already
      curVars = ls(envir=.GlobalEnv)
      if(varName %in% curVars) {
        override = gconfirm(
          sprintf("A variable %s already exists. Overwrite?",varName)
          )
        if(override == FALSE)
          return(TRUE)
      }
      assign_global(make.names(varName),out)
      dispose(.$window) ## clean up
    }
  }
  GUI$makeBody = function(.,container) {
    .$container = container             # store
    ## dispatch various functions depending on type of filename
    if(is.null(.$ext)) {
      tmp = unlist(strsplit(basename(.$filename), split="\\."))
      .$ext = tmp[length(tmp)]
    }
    ## now what is the ext
    switch(.$ext,
           "csv" = .$read_text(sep=","),
           "txt" = .$read_text(sep=""),
           "fwf" = .$read_fwf(sep=","),
           "arff" = .$read_foreign(type="arff"),
           "dbf"= .$read_foreign(type="dbf"),
           "DIF" = .$read_DIF(),
           "dta"= .$read_foreign(type="dta"),
#           "epi"= .$read_foreign(type="epi"),
           "mtp"= .$read_foreign(type="mtp"),
           "octave"= .$read_foreign(type="octave"),
           "sav"= .$read_foreign(type="spss"),
           "ssd"= .$read_foreign(type="ssd"),
           "xport"= .$read_foreign(type="xport"),
           "systat"= .$read_foreign(type="systat"),
#           "xls"= .$read_spreadsheet(type="xls"),
#           "odt" = .$read_spreadsheet(type="odt"),
#           "gnumeric" = .$read_spreadsheet(type="gnumeric"),
           .$read_text(sep=""))         # default
  }
  ## each of these has FUN="string", args=list(), varName
  ## will do do.call(FUN,lapply(args,svalue)) to get answer

  ## ITS ONE OF THESE?
  GUI$read_text = function(.,sep) {
    .$FUN = "read.table"
    .$args = list(file = gedit(.$filename))
    .$allSeps = c(",","\\t","",";","\\s") ## others?
    
    ## see ?read.table for numerous arguments

    g = ggroup(horizontal=FALSE, container=.$container)
    glabel(sprintf("Read %s",basename(.$filename)), container=g)

    tbl <- glayout(container=g)
    tbl[1,1] <- .$AssignToText
    tbl[1,2] <- (.$varName <- gedit("X", container=tbl))
    .$varName[] <- ls(envir=.GlobalEnv)
    visible(tbl) <- TRUE

         
    f= gframe(gettext("Import"), container=g)
    tbl <- glayout(container=f)
    tbl[1,1] <- gettext("header")
    tbl[1,2] <- (.$args[['header']] <- gdroplist(c(TRUE,FALSE), container=tbl))
    tbl[1,3] <- gettext("Skip lines")
    tbl[1,4] <- (.$args[["skip"]] <- gspinbutton(0,1000, container=tbl))
    tbl[2,1] <- gettext("Strip whitespace")
    tbl[2,2] <- (.$args[['strip.white']] <- gdroplist(c(TRUE,FALSE), container=tbl))
    tbl[2,3] <- gettext("Skip blank lines")
    tbl[2,4] <- (.$args[['blank.lines.skip']] <- gdroplist(c(FALSE,TRUE), container=tbl))

    visible(tbl) <- TRUE
    f = gframe(gettext("Attributes"), container=g)
    tbl <- glayout(container=f)
    tbl[1,1] <- gettext("Separator")
    tbl[1,2] <- (.$args[['sep']] <- gedit(sep, container=tbl))
#    tbl[1,2] <- (.$args[['sep']] <- gdroplist(.$allSeps, editable=TRUE,container=tbl))
#    svalue(.$args[['sep']]) <- sep
    
    tbl[1,3] <- gettext("quote")
    tbl[1,4] <- (.$args[['quote']] <- gedit("\"", container=tbl))
    tbl[2,1] <- gettext("Decimal point")
    tbl[2,2] <- (.$args[["dec"]] <- gdroplist(c(".",","), container=tbl))
    tbl[2,3] <- gettext("Comment char.")
    tbl[2,4] <- (.$args[['comment.char']] <- gedit("#", container=tbl))
    tbl[3,1] <- gettext("NA string")
    tbl[3,2] <- (.$args[['na.strings']] <- gedit("NA", container=tbl))

    visible(tbl) <- TRUE

    makePreview = function(...) {
      ## read in
      l <- lapply(.$args, svalue)
      l$nrows = 10
      df= try(do.call(.$FUN,l), silent=TRUE)
      print("DEBUG")
      print(df)
      if(!inherits(df,"try-error")) {
        delete(.$og,.$ig)
        .$ig <- ggroup(horizontal=FALSE, container=.$og, expand=TRUE)
        tmp <- gdf(df,container=.$ig) ## get rownames
##         enabled(tmp) <- FALSE ## too faint
      } else {
        cat(gettext("Error occured:"))
        print(df)
      }
    }

    ## do names?
    f = gframe(gettext("preview"), container=g, expand=TRUE)
    .$og = ggroup(container=f, expand=TRUE)
    .$ig = ggroup(container=.$og, expand=TRUE)                # to be deleted
    makePreview()

    ## now add handler
    sapply(.$args, function(i) addHandlerChanged(i,handler = makePreview))
  }         
  GUI$read_fwf = function(.,sep) {
    .$FUN = "read.fwf"
    .$args = list(file = gedit(.$filename))

    g = ggroup(horizontal=FALSE, container=.$container)
    glabel(paste(gettext("Read"),basename(.$filename),collapse=" "), container=g)

    tbl <- glayout(container=g)
    tbl[1,1] <- .$AssignToText
    tbl[1,2] <- (.$varName <- gedit("X", container=tbl))
    .$varName[] <- ls(envir=.GlobalEnv)
    visible(tbl) <- TRUE

         
    f= gframe(gettext("Import"), container=g)
    tbl <- glayout(container=f)
    tbl[1,1] <- gettext("Header")
    tbl[1,2] <- (.$args[['header']] <- gdroplist(c(FALSE,TRUE), container=tbl))
    tbl[1,3] <- gettext("Separator")
    tbl[1,4] <- (.$args[['sep']] <- gedit(sep, container=tbl))
    tbl[2,1] <- gettext("Skip lines")
    tbl[2,2] <- (.$args[["skip"]] <- gspinbutton(0,1000, container=tbl))
    tbl[2,3] <- gettext("Skip blank lines")
    tbl[2,4] <- (.$args[['blank.lines.skip']] <- gdroplist(c(FALSE,TRUE), container=tbl))
    visible(tbl) <- TRUE
    f = gframe(gettext("Attributes"), container=g)
    tbl <- glayout(container=f)
#    tbl[1,3] <- "quote"
#    tbl[1,4] <- (.$args[['quote']] <- gedit("\"", container=tbl))
#    tbl[2,1] <- "Decimal point"
#    tbl[2,2] <- (.$args[["dec"]] <- gdroplist(c(".",","), container=tbl))
    tbl[1,1] <- gettext("Comment char.")
    tbl[1,2] <- (.$args[['comment.char']] <- gedit("#", container=tbl))
#    tbl[3,1] <- "NA string"
#    tbl[3,2] <- (.$args[['na.strings']] <- gedit("NA", container=tbl))

    visible(tbl) <- TRUE

    ## widths is key here

    f = gframe(gettext("Field widths"), container=g)
    tbl <- glayout(container=f)
    tbl[1,1] <- gettext("widths")
    tbl[1,2] <- (.$args[["widths"]] <- gedit(paste("c(",nchar(readLines(.$filename,n=1)),")",collapse=""), coerce.with=svalue,container=tbl))
    visible(tbl) <- TRUE


    makePreview = function(...) {
      ## read in
      l <- lapply(.$args, svalue)
      l$nrows = 10
      df= try(do.call(.$FUN,l), silent=TRUE)
      if(!inherits(df,"try-error")) {
        delete(.$og,.$ig)
        .$ig <- ggroup(horizontal=FALSE, container=.$og, expand=TRUE)
        tmp <- gdf(df,container=.$ig) ## get rownames
##         enabled(tmp) <- FALSE ## too faint
      } else {
        cat(gettext("Error:"),df,"\n")
      }
    }

    ## do names?
    f = gframe(gettext("preview"), container=g,expand=TRUE)
    .$og = ggroup(container=f, expand=TRUE)
    .$ig = ggroup(container=.$og, expand=TRUE)                # to be deleted
    makePreview()

    ## now add handler
    sapply(.$args, function(i) addHandlerChanged(i,handler = makePreview))

  }
  GUI$read_DIF = function(.) {
    .$FUN = "read.DIF"
    .$args = list(file = gedit(.$filename))

    g = ggroup(horizontal=FALSE, container=.$container)
    glabel(paste(gettext("Read"),basename(.$filename),collapse=" "), container=g)

    tbl <- glayout(container=g)
    tbl[1,1] <- .$AssignToText
    tbl[1,2] <- (.$varName <- gedit("X", container=tbl))
    .$varName[] <- ls(envir=.GlobalEnv)
    visible(tbl) <- TRUE

         
    f= gframe(gettext("Import"), container=g)
    tbl <- glayout(container=f)
    tbl[1,1] <- gettext("Header")
    tbl[1,2] <- (.$args[['header']] <- gdroplist(c(FALSE,TRUE), container=tbl))
    tbl[2,1] <- gettext("Skip lines")
    tbl[2,2] <- (.$args[["skip"]] <- gspinbutton(0,1000, container=tbl))
    tbl[2,3] <- gettext("Skip blank lines")
    tbl[2,4] <- (.$args[['blank.lines.skip']] <- gdroplist(c(FALSE,TRUE), container=tbl))
    tbl[3,1] <- gettext("NA string")
    tbl[3,2] <- (.$args[['na.strings']] <- gedit("NA", container=tbl))
    tbl[3,3] <- gettext("Strings as factors")
    tbl[3,4] <- (.$args[['stringsAsFactors']] <- gdroplist(c(TRUE,FALSE),container=tbl))
    visible(tbl) <- TRUE


    makePreview = function(...) {
      ## read in
      l <- lapply(.$args, svalue)
      l$nrows = 10
      df= try(do.call(.$FUN,l), silent=TRUE)
      if(!inherits(df,"try-error")) {
        delete(.$og,.$ig)
        .$ig <- ggroup(horizontal=FALSE, container=.$og, expand=TRUE)
        tmp <- gdf(df,container=.$ig) ## get rownames
##         enabled(tmp) <- FALSE ## too faint
      } else {
        cat(gettext("Error:"),df,"\n")
      }
    }

    ## do names?
    f = gframe(gettext("preview"), container=g, expand=TRUE)
    .$og = ggroup(container=f, expand=TRUE)
    .$ig = ggroup(container=.$og, expand=TRUE)                # to be deleted
    makePreview()

    ## now add handler
    sapply(.$args, function(i) addHandlerChanged(i,handler = makePreview))

  }    

  GUI$read_foreign = function(.,type) {
    .$FUN = paste("read.",type,sep="",collapse="")
    .$args = list(file=gedit(.$filename)) # all have file as first arg


    fileType = names(.fileExtensions)[sapply(.fileExtensions,function(i) .$ext %in% i)]
        ## strip s
    g = ggroup(horizontal=FALSE, container=.$container)

    glabel(paste(gettext("Read"),basename(.$filename),gettext("as"),popchar(fileType),collapse=" "),
           container=g)
    tbl = glayout(container=g)
    tbl[1,1] <- .$AssignToText
    tbl[1,2] <- (.$varName <- gedit("X", container=tbl))
    .$varName[] <- ls(envir=.GlobalEnv)

    fmls = formals(get(.$FUN))
    nfmls = names(fmls)
    n <- length(nfmls)
    ## add extra arguments without thinking too much
    if(n > 1) {
      for(i in 2:n) {
        tbl[i,1] <- nfmls[i]
        tbl[i,2] <- (.$args[[nfmls[i]]] <-
                     gedit(fmls[[i]], container=tbl,
                           coerce.with = paste("as.",class(fmls[[i]]),sep="", collapse="")
                           ))
        
      }
    }
    
    visible(tbl) <- TRUE

  }

  ## show GUI$show()
  GUI$show()
}
